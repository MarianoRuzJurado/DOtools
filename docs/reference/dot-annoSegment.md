# Annotation modifier for plots

Used for segment the plot for further annotations

## Usage

``` r
.annoSegment(
  object = NULL,
  relSideDist = 0.1,
  aesGroup = FALSE,
  aesGroName = NULL,
  annoPos = "top",
  xPosition = NULL,
  yPosition = NULL,
  pCol = NULL,
  segWidth = 1,
  lty = NULL,
  lwd = 2,
  alpha = NULL,
  lineend = "square",
  annoManual = FALSE,
  mArrow = NULL,
  addBranch = FALSE,
  bArrow = NULL,
  branDirection = 1,
  branRelSegLen = 0.3,
  addText = FALSE,
  textCol = NULL,
  textSize = NULL,
  fontfamily = NULL,
  fontface = NULL,
  textLabel = NULL,
  textRot = 45,
  textHVjust = 0.2,
  hjust = 0.1,
  vjust = -0.5,
  myFacetGrou = NULL,
  aes_x = NULL,
  aes_y = NULL,
  segment_gap = 0.9
)
```

## Arguments

- object:

  ggplot list. Default(NULL).

- relSideDist:

  The relative distance ratio to the y axis range. Default(0.1).

- aesGroup:

  Whether use your group column to add rect annotation.
  Default("FALSE").

- aesGroName:

  The mapping column name. Default(NULL).

- annoPos:

  The position for the annotation to be added. Default("top").

- xPosition:

  The x axis coordinate for the segment. Default(NULL).

- yPosition:

  The y axis coordinate for the segment. Default(NULL).

- pCol:

  The segment colors. Default(NULL).

- segWidth:

  The relative segment width. Default(1).

- lty:

  The segment line type. Default(NULL).

- lwd:

  The segment line width. Default(NULL).

- alpha:

  The segment color alpha. Default(NULL).

- lineend:

  The segment line end. Default("square").

- annoManual:

  Whether annotate by yourself by supplying with x and y coordinates.
  Default(FALSE).

- mArrow:

  Whether add segment arrow. Default(FALSE).

- addBranch:

  Whether add segment branch. Default(FALSE).

- bArrow:

  Whether add branch arrow. Default(FALSE).

- branDirection:

  The branch direction. Default(1).

- branRelSegLen:

  The branch relative length to the segment. Default(0.3).

- addText:

  Whether add text label on segment. Default(FALSE).

- textCol:

  The text colors. Default(NULL).

- textSize:

  The text size. Default(NULL).

- fontfamily:

  The text fontfamily. Default(NULL).

- fontface:

  The text fontface. Default(NULL).

- textLabel:

  The text textLabel. Default(NULL).

- textRot:

  The text angle. Default(NULL).

- textHVjust:

  The text distance from the segment. Default(0.2).

- hjust:

  The text hjust. Default(NULL).

- vjust:

  The text vjust. Default(NULL).

- myFacetGrou:

  Your facet group name to be added with annotation when object is a
  faceted object. Default(NULL).

- aes_x:

  = NULL You should supply the plot X mapping name when annotate a
  facetd plot. Default(NULL).

- aes_y:

  = NULL You should supply the plot Y mapping name when annotate a
  facetd plot. Default(NULL).

- segment_gap:

  = 0.9 define the gap between segmentation brackets

## Value

ggplot

## Author

Mariano Ruz Jurado (edited from: Jun Zhang)
