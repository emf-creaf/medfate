# An ongoing project

R package **medfate** and its associated R packages should be viewed as an ongoing research project for the development of forest ecosystem modelling tools. The **medfate** R package is the result of collaborative work between modellers and experts in different disciplines. Since successful modelling projects involve long-term investments and the participation of multiple teams, we are open to further expanding the set of people contributing to the project. Normally, contributors will start participating as model users, but may soon have their own ideas in how to improve the model or encounter some issues to be solved. 

# Contributing to medfate

Contributions to the development of **medfate** can be done in different aspects:

 1. **Model design and formulation**: If your expertise includes any of the processes that are modelled in the package and you feel that your expertise could be helpful to improve the package, you are more than welcome to contact us! If you are familiarized with Git, GitHub and R package development, you can fork the package, make changes and then a pull request (see below). Otherwise, other forms of collaboration can be established. While contributions are welcome, we do not want to have multiple, diverging, versions of the simulation models. Hence, we want to centralize and review modifications, so that former package functionality is not lost.
 2. **Model parameterization**: Finding suitable parameter values for trait-based models is hard, and requires gathering data from multiple databases. Efforts to find species parameter values required for **medfate** can be made available to others by including new species parameter tables, such as `SpParamsMED`. We are currently developing on a companion package called [**medfatetraits**](https://github.com/emf-creaf/medfatetraits/) that should be helpful to define and populate new species parameter tables.
 3. **Model evaluation**:  Simulation models should be tested extensively, and there is a lot to be done in this respect in the case of **medfate**. Hence, we will appreciate help in this area, for example pointing at interesting validation data sets. They should lead to new package vignettes showing the performance of the model in different situations.

## Reporting bugs and suggesting enchancements

If you want to report a bug or suggest an enhancement, it's a good idea to file an issue to the medfate repository at GitHub. If you’ve found a bug, please file an issue that illustrates the bug with a minimal 
[reprex](https://www.tidyverse.org/help/#reprex).

### Code contributions

Before making contributions to the package R or C++ code, make sure someone from the **medfate** team agrees that the change you suggest is needed. 

*   Fork the package and clone onto your computer. If you haven't done this before, we recommend using `usethis::create_from_github("emf-creaf/medfate", fork = TRUE)`.

*   Install all development dependences with `devtools::install_dev_deps()`, and then make sure the package passes R CMD check by running `devtools::check()`. 
    If R CMD check doesn't pass cleanly, it's a good idea to ask for help before continuing. 
*   Create a Git branch for your pull request (PR). We recommend using `usethis::pr_init("brief-description-of-change")`.

*   Make your changes, commit to git, and then create a PR by running `usethis::pr_push()`, and following the prompts in your browser.
    The title of your PR should briefly describe the change.
    The body of your PR should contain `Fixes #issue-number`.

*  For user-facing changes, add a bullet to the top of `NEWS.md` (i.e. just below the first header). Follow the style described in <https://style.tidyverse.org/news.html>.
