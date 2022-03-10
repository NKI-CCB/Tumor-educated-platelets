library(renv)

#From https://rstudio.github.io/renv/articles/renv.html

#Call renv::init() to initialize a new project-local environment with a private R library,
renv::init()

#Work in the project as normal, installing and removing new R packages as they are needed in the project,

#Call renv::snapshot() to save the state of the project library to the lockfile (called renv.lock),
renv::snapshot()

#Continue working on your project, installing and updating R packages as needed.

#Call renv::snapshot() again to save the state of your project library if your attempts to update R packages were successful,
renv::snapshot()

#or call renv::restore() to revert to the previous state as encoded in the lockfile if your attempts to update packages introduced some new problems.
renv::restore()
