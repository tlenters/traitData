# Start server
bookdown::serve_book(dir = ".", output_dir = "_book", preview = TRUE, in_session = TRUE, quiet = FALSE)

# Stop server
servr::daemon_stop(1)
