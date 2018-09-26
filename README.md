# landscape.windows
 Description: These functions create landscape windows into an area defined by a shapefile (shape format). There are two different methods: regular and random. Using regular method: latitudinal transects aligned every n kilometers are used to place the center of each windows. Windows are aligned along these transects at n kilometers distances between their centers. Using random method: n windows are placed randomly in the defined area. With default values, to be validated, a window must contain at least 25 sites distributed over at least 80% of the window’s area (at least 1 site in 8 cells of the 3 x 3 grid).

## Installation

```r
devtools::install_github("tleboucher/landscape.windows")
```
