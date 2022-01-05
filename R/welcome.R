.onAttach <- function(lib, pkg)  {
  packageStartupMessage("This is medfate [ver. ",
                        utils::packageDescription("medfate",
                                                  fields="Version"),"]",
                        appendLF = TRUE)
}
