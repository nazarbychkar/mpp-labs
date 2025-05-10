# === 1. Задаємо числовий ряд === 
series <- c(3.2, 7.8, 1.5, 9.0, 4.6) 

# === 2. Визначення функцій === 

# Функція для сортування числового ряду 
sort_series <- function(series) { 
  sort(series) 
}

# Побудова інтервалів на основі пуассонівського розподілу 
build_poisson_intervals <- function(series, alphabet_size, lambda = 1) {
  min_val <- min(series)
  max_val <- max(series)
  
  # Створюємо вектор ймовірностей за розподілом Пуассона
  x <- 0:(alphabet_size-1)
  poisson_probs <- dpois(x, lambda)
  poisson_probs <- poisson_probs / sum(poisson_probs)  # Нормалізація ймовірностей
  
  # Обчислюємо кумулятивні ймовірності
  cum_probs <- cumsum(poisson_probs)
  
  # Створюємо межі інтервалів на основі кумулятивних ймовірностей
  range <- max_val - min_val
  breaks <- min_val + cum_probs * range
  breaks <- c(min_val, breaks)
  
  return(breaks)
}

# Перетворення чисел у відповідні символи алфавіту
map_to_symbols <- function(series, breaks, alphabet) {
  indices <- findInterval(series, breaks, rightmost.closed = TRUE)
  indices[indices == 0] <- 1
  indices[indices > length(alphabet)] <- length(alphabet)
  return(alphabet[indices])
}

# Побудова матриці передування символів
build_transition_matrix <- function(symbols, alphabet) {
  n <- length(alphabet)
  matrix <- matrix(0, nrow = n, ncol = n, dimnames = list(alphabet, alphabet))
  
  for (i in 1:(length(symbols) - 1)) {
    from <- symbols[i]
    to <- symbols[i + 1]
    matrix[from, to] <- matrix[from, to] + 1
  }
  
  return(matrix)
}

# Головна функція
run_lab <- function(series, alphabet, lambda = 1) {
  alphabet_size <- length(alphabet)
  sorted <- sort_series(series)
  
  # Використовуємо розподіл Пуассона для створення інтервалів
  breaks <- build_poisson_intervals(sorted, alphabet_size, lambda)
  
  symbols <- map_to_symbols(series, breaks, alphabet)
  transition_matrix <- build_transition_matrix(symbols, alphabet)
  
  list(
    series = series,
    sorted = sorted,
    breaks = breaks,
    symbols = symbols, 
    matrix = transition_matrix
  )
}

# === 3. Визначення алфавіту === 
alphabet <- c("A", "B", "C", "D", "E")

# === 4. Запуск алгоритму з вимірюванням часу виконання === 
start_time <- Sys.time()
result <- run_lab(series, alphabet, lambda = 1.5)  # Використовуємо параметр lambda = 1.5 для прикладу
end_time <- Sys.time()
execution_time <- end_time - start_time

# === 5. Вивід результатів === 
cat("\nВхідний числовий ряд:\n")
cat(result$series, sep = ", ")

cat("\n\nВідсортований числовий ряд:\n")
cat(result$sorted, sep = ", ")

cat("\n\nМежі інтервалів за розподілом Пуассона:\n")
cat(result$breaks, sep = ", ")

cat("\n\nЛінгвістичний ряд:\n")
cat(result$symbols, sep = " ")

cat("\n\nМатриця передування:\n")
print(result$matrix)

cat("\nЧас виконання (мілісекунди):\n")
cat(round(as.numeric(execution_time, units = "secs") * 1000, 3))
