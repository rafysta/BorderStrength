packages <- c("optparse", "data.table")

for (pkg in packages) {
  if (!require(pkg, character.only = TRUE)) {
    install.packages(pkg, repos = "https://cloud.r-project.org", dependencies = TRUE)
  }
}
