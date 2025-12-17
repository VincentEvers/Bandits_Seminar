
library(dplyr)
library(tidyr)
library(ggplot2)


set.seed(42)

d <- 3
T <- 500

eta <- sqrt(log(d)/(T*d))

w <- rep(1/d, d)

cost_arm_1 <- pmax(pmin(rnorm(T, mean = 0.3, sd = 0.2), 1), 0)
cost_arm_2 <- pmax(pmin(rnorm(T, mean = 0.5, sd = 0.2), 1), 0)
cost_arm_3 <- pmax(pmin(rnorm(T, mean = 0.7, sd = 0.2), 1), 0)
cost_mat <- matrix(data = c(cost_arm_1, cost_arm_2, cost_arm_3), nrow = T, ncol = d)

w_history <- matrix(NA, nrow = T, ncol = d)
arm_chosen <- numeric(T)

sample_arm <- function(weights) {
  sample(1:d, size = 1, prob = weights)
}

for (t in 1:T) {
  pt <- sample_arm(w)
  
  loss <- cost_mat[t, pt]
  
  w_tilde <- w
  w_tilde[pt] <- w[pt] * exp(-eta * loss / w[pt])
  w <- w_tilde / sum(w_tilde)
  
  arm_chosen[t] <- pt
  
  w_history[t, ] <- w
}

# Weight Plot

w_df <- as.data.frame(w_history)
w_df$round <- 1:nrow(w_df)

w_long <- w_df %>%
  pivot_longer(cols = -round,
               names_to = "arm",
               values_to = "weight") %>%
  mutate(arm = factor(arm))

sim_plot_stacked <- ggplot(w_long, aes(x = round, y = weight, fill = arm)) +
  geom_area(alpha = 0.8, linewidth = 0.25, color = "white") +
  scale_fill_manual(values = c("#C71F37", "#F26A1B", "#F7C12A"),
                    labels = c("1 ~ N(0.3, 0.04)", "2 ~ N(0.5, 0.04)", "3 ~ N(0.7, 0.04)")) +
  labs(
    title = "Evolution of weight vector w",
    x = "Round",
    y = "Weight",
    fill = "Arm"
  ) +
  theme_minimal(base_size = 14) +
  theme(plot.title = element_text(hjust = 0.5))

ggsave("Code_Sim_Weight.png", sim_plot_stacked, width = 12, height = 6, dpi = 800)


## Line Plot


cost_df <- as.data.frame(cost_mat)
cost_df <- cbind(cost_df, cost_mat[cbind(1:T, arm_chosen)])

cost_df <- cost_df - cost_mat[, 1]

cost_df <- as.data.frame(apply(cost_df, 2, cumsum))
cost_df <- cbind(1:T, cost_df)
cost_df <- cbind(cost_df, 2 * sqrt(d*log(d)*1:T))
colnames(cost_df) <- c("Round", "Arm 1", "Arm 2", "Arm 3", "Our algorithm", "Regret bound")

cost_df <- pivot_longer(cost_df, cols = -Round, names_to = "Arm", values_to = "Regret")

sim_plot_line <- ggplot(cost_df, aes(x = Round, y = Regret, color = Arm)) +
  geom_line(linewidth = 1.5, alpha = 0.8) +
  scale_color_manual(values = c("#C71F37", "#F26A1B", "#F7C12A", "#457B9D", "#6A4C93")) +
  labs(title = "Regret for every arm and our strategy") +
  theme_minimal(base_size = 14) +
  theme(plot.title = element_text(hjust = 0.5))

ggsave("Code_Sim_Regret.png", sim_plot_line, width = 12, height = 6, dpi = 800)


