import pandas as pd
import matplotlib.pyplot as plt
import re
from typing import Union
import os
def format_terms_gsea(df, term_col, cutoff=35):
    def remove_whitespace_around_newlines(text):
        # Replace whitespace before and after newlines with just the newline
        return re.sub(r'\s*\n\s*', '\n', text)

    newterms = []
    for text in df[term_col]:
        newterm, text_list_nchar, nchar, limit = [], [], 0, cutoff
        text_list = text.split(' ')
        for txt in text_list:  # From text_list get a list where we sum nchar from a word + previous word
            nchar += len(txt)
            text_list_nchar.append(nchar)
        for idx, word in enumerate(text_list_nchar):
            if word > limit:  # If we have more than cutoff characters in len add a break line
                newterm.append('\n')
                limit += cutoff
            newterm.append(text_list[idx])
        newterm = ' '.join(newterm)
        cleanterm = remove_whitespace_around_newlines(newterm)  # remove whitespace inserted
        newterms.append(cleanterm)
    df[term_col] = newterms

    return df


def split_bar_gsea(df: pd.DataFrame,
                   term_col: str,
                   col_split: str,
                   cond_col: str,
                   pos_cond: str,
                   cutoff: int = 40,
                   log10_transform: bool = True,
                   figsize: tuple = (12, 8),
                   topN: float = 10,
                   colors_pairs: list = ['sandybrown', 'royalblue'],
                   alpha_colors: float = 0.3,
                   path: Union[None, str] = None,
                   spacing: float = 5,
                   txt_size: float = 12,
                   filename: str = 'SplitBar.svg',
                   title: str = 'Top 10 GO Terms in each Condition',
                   show: bool = True) -> Union[None, plt.axis]:
    """ **Split BarPlot for GO terms**

    This function generates a split barplot. This is a plot where the top 10 Go terms
    are shown, sorted based on a column ('col_split'). Two conditions are shown at the same
    time. One condition is shown in the positive axis, while the other in the negative one.
    The condition to be shown as positive is set with 'pos_col'.

    The GO terms will be shown inside the bars, if the term is too long, using 'cutoff',
    you can control the maximum number of characters per line.

    !! Pre-filter the dataframe to contain significant Terms

    :param df: dataframe with the results of a gene set enrichment analysis
    :param term_col: column in the dataframe that contains the terms
    :param col_split: column in the dataframe that will be used to sort and split the plot
    :param cond_col: column in the dataframe that contains the condition information
    :param pos_cond: condition that will be shown in the positive side of the plot
    :param cutoff: maximum number of characters per line
    :param log10_transform: if col_split contains values between 0 and 1, assume they are pvals and apply a -log10 transformation
    :param figsize: figure size
    :param topN: how many terms are shown
    :param path: path to save the plot
    :param filename: filename for the plot
    :param spacing: space to add between bars and origin. It is a percentage value, indicating that the bars start at 5 % of the maximum X axis value.
    :param txt_size: size of the go terms text
    :param alpha_colors: alpha value for the colors of the bars
    :param colors_pairs: colors for each condition (1st color --> negative axis; 2nd color --> positive axis)
    :param title: title of the plot
    :param show: if False, the axis is return
    :return: None or the axis
    """

    if len(df[cond_col].unique()) != 2:
        if len(df[cond_col].unique()) > 2:
            assert len(df[cond_col].unique()) == 2, 'Not implement - Only 1 or 2 conditions can be used'
        elif len(df[cond_col].unique()) == 1:
            print('!!! WARNING - There are no terms for one of the conditions')
        else:
            assert len(df[cond_col].unique()) == 2, 'Not implement - Only 1 or 2 conditions can be used'

    print('!!! Assuming GO Terms are preprocessed (Only Significant terms included)')

    df = df.copy()  # Ensure we do not modify the input
    jdx = list(df.columns).index(cond_col)  # Get index of the condition column

    # Update the col_split values; Positive values for one condition and
    # negative for the other positive. The positive is set by the 'pos_cond' argument
    min_val, max_val = df[col_split].min(), df[col_split].max()
    is_pval = True if (min_val >= 0) and (max_val <= 1) else False
    if is_pval and log10_transform:
        print('Assuming col_split contains Pvals, apply -log10 transformation')
        df['-log10(Padj)'] = -np.log10(df[col_split])
        col_split = '-log10(Padj)'
        spacing = .5  # Correct spacing in case it was not specified
    df[col_split] = [val if df.iloc[idx, jdx] == pos_cond else -val for idx, val in
                     enumerate(df[col_split])]  # Set negative and positive values for each condition

    # Format the Terms
    df[term_col] = df[term_col].str.capitalize()  # Capitalise
    df = format_terms_gsea(df, term_col, cutoff)  # Split terms too long in several rows

    # Get the dataframe for the positive and negative axis
    df_pos = df[df[cond_col] == pos_cond].sort_values(col_split, ascending=False).head(int(topN))
    df_neg = df[df[cond_col] != pos_cond].sort_values(col_split).head(int(topN))

    # Check that the size of the dataframes is equal
    if len(df_pos) != len(df_neg):
        print('Different number of GO Terms in positive and negative axis, adding empty rows')
        print(f'Positive side has {len(df_pos)} and Negative side has {len(df_neg)}')
        missing_rows = topN - len(df_pos) if len(df_pos) < len(df_neg) else topN - len(df_neg)
        missing_rows_data = [np.nan for val in range(len(df_pos.columns))]
        missing_df = pd.DataFrame([missing_rows_data] * missing_rows, columns=list(df_pos.columns))
        missing_df[term_col] = ''
        missing_df[col_split] = 0
        if len(df_pos) > len(df_neg):
            df_neg = pd.concat([df_neg, missing_df])
        else:
            df_pos = pd.concat([df_pos, missing_df])

    spacing_unit = np.abs(df[col_split]).max() * spacing / 100
    # Actual Plot
    fig, axs = plt.subplots(1, 1, figsize=figsize)
    y_pos = range(int(topN))

    # Plot bars for "Down" condition (positive values) on the left side
    bars_down = axs.barh(y_pos,
                         df_neg[col_split].sort_values(ascending=False),
                         left=-spacing_unit, color=colors_pairs[0],
                         align='center', alpha=alpha_colors)

    # Plot bars for "Up" condition (negative values) on the right side
    bars_up = axs.barh(y_pos,
                       df_pos[col_split].sort_values(),
                       left=spacing_unit, color=colors_pairs[1],
                       align='center', alpha=alpha_colors)

    # Layout
    axs.spines[['left', 'top', 'right']].set_visible(False)
    axs.set_yticks([])
    axs.set_xlim(-np.abs(df[col_split]).max(), np.abs(df[col_split]).max())
    axs.set_xlabel(col_split, fontsize=18)
    axs.set_title(title, fontsize=20)
    axs.grid(False)
    plt.vlines(0, -1, float(topN) - .5, color='k', lw=1)
    axs.set_ylim(-.5, float(topN))

    # Add text labels for each bar (GO term name)
    for i, bar in enumerate(bars_up):
        # Add the GO term for "Up" bars (positive)
        axs.text(spacing_unit * 2, bar.get_y() + bar.get_height() / 2,
                 df_pos.sort_values(col_split)[term_col].iloc[i],
                 va='center', ha='left', color='k', fontweight='bold', fontsize=txt_size)

    for i, bar in enumerate(bars_down):
        # Add the GO term for "Down" bars (negative)
        axs.text(-spacing_unit * 2, bar.get_y() + bar.get_height() / 2,
                 df_neg.sort_values(col_split, ascending=False)[term_col].iloc[i],
                 va='center', ha='right', color='k', fontweight='bold', fontsize=txt_size)
    # Save Plot
    if path is not None:
        plt.savefig(os.path.join(path, filename), bbox_inches='tight')

    # If show is False, return axs
    if not show:
        return axs
    else:
        return

