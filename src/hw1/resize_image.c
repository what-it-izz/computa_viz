#include <math.h>
#include "image.h"

// performs nearest neighbor interpolation
float nn_interpolate(image im, int c, float h, float w)
{
    return get_pixel(im, c, round(h), round(w));
}

image nn_resize(image im, int h, int w)
{
    // Create a new image that is h x w and the same number of channels as im
    image nn_resized = make_image(im.c,h,w);
    
    float w_ratio, h_ratio, old_w, a, b, old_h, curr_pixel;
    h_ratio = (float) im.h / (float) h;
    w_ratio = (float) im.w / (float) w;
    a = h_ratio * 0.5 - 0.5;
    b = w_ratio * 0.5 - 0.5;

    // Loop over the pixels and map back to the old coordinates
    // Use nearest-neighbor interpolate to fill in the image
    for (int c = 0; c < im.c; c++) {
        for (int i = 0; i < h; i++) {
            old_h = i * h_ratio + a;
            for (int j = 0; j < w; j++) {
                old_w = j * w_ratio + b;
                curr_pixel = nn_interpolate(im, c, old_h, old_w);
                set_pixel(nn_resized, c, i, j, curr_pixel);
            }
        }
    }
    return nn_resized;
}

// performs bilinear interpolation
float bilinear_interpolate(image im, int c, float h, float w)
{
    float floor_h, ceil_h, floor_w, ceil_w, v1, v2, v3, v4, d1, d2, d3, d4, q1, q2, q;
    floor_h = floor(h);
    floor_w = floor(w);
    ceil_h = ceil(h);
    ceil_w = ceil(w);
    
    v1 = get_pixel(im, c, ceil_h, floor_w);
    v2 = get_pixel(im, c, ceil_h, ceil_w);
    v3 = get_pixel(im, c, floor_h, floor_w);
    v4 = get_pixel(im, c, floor_h, ceil_w);
    d1 = w - floor_w;
    d2 = ceil_w - w;
    d3 = ceil_h - h;
    d4 = h - floor_h;

    q1 = v1*d2 + v2*d1;
    q2 = v3*d2 + v4*d1;
    q = q1*d4 + q2*d3;

    return q;
}

image bilinear_resize(image im, int h, int w)
{
    // Create a new image that is h x w and the same number of channels as im
    image bl_resized = make_image(im.c,h,w);
    
    float w_ratio, h_ratio, old_w, a, b, old_h, curr_pixel;
    h_ratio = (float) im.h / (float) h;
    w_ratio = (float) im.w / (float) w;
    a = h_ratio * 0.5 - 0.5;
    b = w_ratio * 0.5 - 0.5;

    // Loop over the pixels and map back to the old coordinates
    // Use bilinear interpolate to fill in the image
    for (int c = 0; c < im.c; c++) {
        for (int i = 0; i < h; i++) {
            old_h = i * h_ratio + a;
            for (int j = 0; j < w; j++) {
                old_w = j * w_ratio + b;
                curr_pixel = bilinear_interpolate(im, c, old_h, old_w);
                set_pixel(bl_resized, c, i, j, curr_pixel);
            }
        }
    }
    return bl_resized;
}

