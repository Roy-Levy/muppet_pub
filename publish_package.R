# publish_package_safe.R
publish_package <- function(pkg_path = ".", docs_path = "docs") {

  if (!dir.exists(pkg_path)) stop("Package path not found!")

  old_wd <- setwd(pkg_path)
  on.exit(setwd(old_wd), add = TRUE)

  # Build pkgdown site if pkgdown is installed
  if (requireNamespace("pkgdown", quietly = TRUE)) {
    message("Building pkgdown site...")
    pkgdown::build_site(pkg_path)
  } else {
    warning("pkgdown not installed. Skipping site build.")
  }

  # Check if there are any changes to commit
  status <- system("git status --porcelain", intern = TRUE)
  if (length(status) == 0) {
    message("No changes to commit. Exiting.")
    return(invisible(NULL))
  }

  # Tag the previous commit as a backup
  old_commit <- system("git rev-parse HEAD", intern = TRUE)
  tag_name <- paste0("backup-", substr(old_commit, 1, 7))
  system(paste("git tag -f", tag_name))
  message(paste("Tagged previous commit as", tag_name))

  # Ask for a commit message
  commit_msg <- readline("Enter commit message: ")

  # Stage and commit changes
  system("git add .")
  system(paste("git commit -m", shQuote(commit_msg)))

  # Force push to GitHub
  system("git push --force")

  message("✅ Package and website updated on GitHub!")
  message(paste("Previous commit tagged as:", tag_name))
}

# Run immediately
publish_package()
