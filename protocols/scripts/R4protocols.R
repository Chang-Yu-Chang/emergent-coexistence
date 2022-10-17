

#' Draw the 96-well plate layout
#' @description Draw the 96 well plate layout for experiment
#' @name plate_layout_draw
#' @param well_label label for each well.
#' @param row_position numeric vector. the desired row position.
#' @param col_position numeric vector. the desired column position.
#' @param annotation logical. If FALSE there is no annotation on the plate.
#' @param annotate_each_well logical. If FALSE (the default) wells with the same attributes will share one annotated text
#' in the center of all such wells, otherwise each well has its own annotated text.
#' @param well_size desired size of wells.
#' @param text_size desired size of annotated texts.
#' @param byrow logical. If FALSE (the default) the matrix is filled by columns.
#' @return ggplot object
#' @examples
#'  plate_layout_draw(
#'  label = factor(rep(c(1:6), 3), 1:6),
#'  row_position = rep(1:6, 3),
#'  col_position = rep(1:3, each = 6)
#'  )
#' \dontrun{
#' hello("test")
#' }
#' @export

plate_layout_draw <- function(
  label = c("C", "E", "P", "R", "rCEPR"),
  row_position = rep(2, 5),
  col_position = seq(2, 10, 2),
  annotation = FALSE,
  annotate_each_well = FALSE,
  well_size = 5,
  text_size = 2,
  byrow = FALSE,
  ...
) {

  # Check ----
  stopifnot(length(label) == length(row_position) |
      length(label) == length(col_position))

  # Format input label ----
  if (is.numeric(label)) {
    plate <-
      data.frame(
        Label = label,
        Row = row_position,
        Col = col_position
      )
  } else if (is.character(label)) {
    plate <-
      data.frame(
        Label = factor(label, level = unique(label)),
        Row = row_position,
        Col = col_position
      )
  } else if (is.factor(label)) {
    plate <-
      data.frame(
        Label = label,
        Row = row_position,
        Col = col_position
      )
  }

  # Empty well ----
  empty_well <- data.frame(
    Row = rep(1:8, each = 12),
    Col = rep(1:12, 8)
  )

  # Plot blank plate
  if (is.null(label)) {
    return(
    ggplot() +
      geom_point(data = empty_well, aes(x = Col, y = Row), shape = 1, size = well_size, color = "black", fill = NA) +
      scale_x_continuous(name = "", breaks = 1:12, labels = 1:12, limits = c(1, 12), position = "top") +
      scale_y_reverse(name = "", lim = c(8, 1), breaks = 1:8, labels = LETTERS[1:8]) +
      theme_bw() +
      NULL
    )
  }

  # Plot segments
  g_empty_well <- geom_point(data = empty_well, aes(x = Col, y = Row),
    shape = 1, size = well_size, color = "black", fill = NA)

  g_plate <- geom_point(data = plate, aes(x = Col, y = Row, fill = Label),
    pch = 21, size = well_size, show.legend = FALSE)


  # Plot ----
  if (annotation == TRUE) {
    # Annotation
    ann_pool <- annotate("text", x = tapply(plate$Col, plate$Label, mean), y = tapply(plate$Row, plate$Label, mean),
      label = unique(plate$Label), size = text_size)
    ann_indi <- annotate("text", x = plate$Col, y = plate$Row, label = plate$Label, size = text_size)

    if (annotate_each_well == TRUE) {
      ggplot() + g_plate + g_empty_well + ann_indi +
        scale_x_continuous(name = "", breaks = 1:12, labels = 1:12, limits = c(1, 12), position = "top") +
        scale_y_reverse(name = "", lim = c(8, 1), breaks = 1:8, labels = LETTERS[1:8]) +
        theme_bw() +
        NULL

    } else if (annotate_each_well == FALSE) {
      ggplot() + g_plate + g_empty_well + ann_pool +
        scale_x_continuous(name = "", breaks = 1:12, labels = 1:12, limits = c(1, 12), position = "top") +
        scale_y_reverse(name = "", lim = c(8, 1), breaks = 1:8, labels = LETTERS[1:8]) +
        theme_bw() +
        NULL

    }
  } else if (annotation == FALSE) {
    # Only plate
    ggplot() + g_plate + g_empty_well +
      scale_x_continuous(name = "", breaks = 1:12, labels = 1:12, limits = c(1, 12), position = "top") +
      scale_y_reverse(name = "", lim = c(8, 1), breaks = 1:8, labels = LETTERS[1:8]) +
      theme_bw() +
      NULL
  }
}


