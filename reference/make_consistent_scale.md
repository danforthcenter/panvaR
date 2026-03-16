# Generate discrete scale

Make a ggplot discrete scale that has consistent correspondence between
aesthetics and variables

## Usage

``` r
make_consistent_scale(
  values,
  vars,
  type = c("fill", "color", "shape"),
  show.example = F,
  name = NULL
)
```

## Arguments

- values:

  vector, values to supply to given scale (e.g. colors or shapes)

- vars:

  character, variables to link to values. Will have one to one
  correspondence to values based on order

- type:

  character, either "fill", "color" or "shape"

- show.example:

  boolean, print a plot to show what values are tied to which variables

- name:

  optional, string to name scale, will appear as legend heading when
  used in a ggplot

## Value

ggplot scale object for use in downstream plotting

## Examples

``` r
# work in progress

# fill
my.scale <- 
make_consistent_scale(
values =  c("red", "blue", "gold"),
vars = c("A", "B", "C"),
show.example = TRUE)


# shape 
my.scale <- 
make_consistent_scale(
values =  c(21, 22, 23),
vars = c("A", "B", "C"),
type = "shape",
show.example = TRUE)
```
