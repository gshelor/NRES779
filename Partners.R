
generate_partner_schedule <- function(names, week_number, past_assignments = NULL) {
  set.seed(2024)
  # Check if the number of names is even or odd
  num_names <- length(names)
  num_pairs <- num_names %/% 2
  has_odd_student <- num_names %% 2 != 0
  
  # Initialize an empty matrix to store partner assignments
  partner_matrix <- matrix(NA, nrow = num_pairs, ncol = 2)
  
  # Create a vector to keep track of which students have been assigned partners in the current week
  assigned <- rep(FALSE, length(names))
  
  # If there are past assignments, mark the students who have already been assigned partners
  if (!is.null(past_assignments)) {
    for (pair in past_assignments) {
      assigned[pair] <- TRUE
    }
  }
  
  # Generate random partner assignments until all pairs have been formed at least once
  while (anyNA(partner_matrix)) {
    # Shuffle the names
    shuffled_names <- sample(names)
    
    # Assign partners for the current week
    for (i in 1:num_pairs) {
      partner_matrix[i, ] <- shuffled_names[c(i, i + num_pairs)]
      assigned[which(names == partner_matrix[i, 1])] <- TRUE
      assigned[which(names == partner_matrix[i, 2])] <- TRUE
    }
    
    # Check if any pair has been assigned twice or if any student has been assigned partners in consecutive weeks
    if (!any(duplicated(partner_matrix)) && !any(assigned)) {
      break
    }
  }
  
  # If there's an odd number of students, add one group of three
  if (has_odd_student) {
    # Find the odd student
    odd_student <- setdiff(names, unlist(partner_matrix))
    
    # Randomly select a pair to form a group of three
    idx <- sample(1:num_pairs, 1)
    partner_matrix <- rbind(partner_matrix, c(partner_matrix[idx, 2], odd_student))
    partner_matrix[idx, 2] <- odd_student
  }
  
  # Return the partner assignments for the given week number
  partner_matrix
}

# Example usage:
names=c("jared", "colton",
        "martin", "kushi",
        "hannah_g", "elly",
        "griffin_p",
        "thomas", "faith",
        "griffin_s","kadin",
        "corina", "victoria",
        "alisa", "rachel",
        "sage", "hannah_p",
        "jordan","anita")

# Initialize an empty list to store partner schedules for each week
partner_schedules <- list()

# Generate partner schedules for 5 weeks
for (week_number in 1:5) {
  if (week_number == 1) {
    past_assignments <- NULL
  } else {
    past_assignments <- cbind(partner_schedules[[week_number - 1]])
  }
  partner_schedule <- generate_partner_schedule(names, week_number, past_assignments)
  partner_schedules[[week_number]] <- partner_schedule
}

# Print partner schedules for each week
for (week_number in 1:5) {
  cat("Week", week_number, ":\n")
  print(partner_schedules[[week_number]])
  cat("\n")
}

