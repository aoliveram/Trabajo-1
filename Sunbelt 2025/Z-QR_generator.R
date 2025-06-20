#install.packages("qrcode")

library(qrcode)
qr_code("https://github.com/aoliveram/Trabajo-1/tree/main/Sunbelt%202025") |>
  generate_svg("diffusion_of_innovations_heatmaps.svg", background = "transparent")
