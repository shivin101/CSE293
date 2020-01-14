function img= linearize_val(img)

img = img-min(min(img));
img = img/max(max(img));
img = 255*(img);

end