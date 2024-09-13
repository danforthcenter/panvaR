To have `PanvaR` in your dev environment please use either:-

1. [`remotes`](https://cran.r-project.org/web/packages/remotes/index.html)

2. [`devtools`](https://cran.r-project.org/web/packages/devtools/index.html)

The above package are the responsibility of the end user.

After you have either `devtools` or `remotes` installed use either:-

With `remotes` :- `remotes::install_github("VectorFrankenstein/panvaR")`

or

With `devtools` :- `devtools::install_github("VectorFrankenstein/panvaR")`

After this you should be ready to play around with `PanvaR`.

Having a fresh `Conda`/`mamba`/`Renv` could make from a smoother experience install/using `PanvaR`, but `PanvaR` does have a fairly large dependency profile on the `tidyverse` and being close to the latest version of `tidyverse` will make for a smooth user experience.