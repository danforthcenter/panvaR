To have `PanvaR` in your dev environment please use the following steps:-

1. Install [`remotes`](https://cran.r-project.org/web/packages/remotes/index.html)

2. After you have `remotes` installed use:- `remotes::remotes::install_github("danforthcenter/panvaR", build_vignettes = TRUE)`

After this you should be ready to play around with `PanvaR`.

Having a fresh `Conda`/`mamba`/`Renv` could make from a smoother experience install/using `PanvaR`. `PanvaR` does have a fairly large dependency profile on the `tidyverse` and being close to the latest version of `tidyverse` will make for a smooth user experience.

Once you have `panvaR` installed you can use `browseVignettes("panvaR")` to get started.
