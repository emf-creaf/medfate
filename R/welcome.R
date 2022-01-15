.onAttach <- function(lib, pkg)  {
  packageStartupMessage("Package 'medfate' [ver. ",
                        utils::packageDescription("medfate",
                                                  fields="Version"),"]",
                        appendLF = TRUE)
}
