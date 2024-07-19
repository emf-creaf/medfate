library(hexSticker)

logo_image <- fs::path("data-raw", "medfate_logo_image.png")
sticker(
  logo_image,
  package = "medfate", p_size = 18, p_y = 1.64, p_color = "#D4AA00",
  s_x = 0.93, s_y = .85, s_width = .58,
  filename = fs::path("data-raw", "medfate.png"),
  #   url = "emf.creaf.cat", u_size = 6, u_color = "#BFD77A", u_y = .2, u_x = 1.2,
  h_fill = "#6C5353", h_color = "#D4AA00"
)
