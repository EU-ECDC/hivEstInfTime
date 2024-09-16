# 1. Update local packages -------------------------------------------------------------------------
dir.create('library', showWarnings = FALSE)
.libPaths('./library')
pak::local_install_deps(root = '.', dependencies = 'hard', upgrade = TRUE)
