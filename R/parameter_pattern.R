# The below two functions are used to define regular expression

# This function is used for a variable name, avoiding partial matches (e.g., "CE" does not match "CEN")
# The variable must appear immediately after "." or "_"
# and must be followed by ".", "$" (for thresholds), or nothing
parameter_pattern <- function(parameter) {
  paste0(
    "(\\.|_)",
    parameter,
    "(\\.|\\$|$)"
  )
}

# This function is used for a complete string such as "F1.by.X1"
# The string must appear immediately after "_" and nothing follows it
parameter_string_pattern <- function(parameter_string) {
  parameter_string <- gsub("\\.", "\\\\.", parameter_string)
  paste0("_", parameter_string, "$")
}
