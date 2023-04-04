#include <stdio.h>
#include <string.h>
#include <assert.h>
#include <math.h>
#include "image.h"

/* HELPER FUNCS */
// returns a pointer to the pixel at the correct index of an image given c, h, & w
float* get_pixel_helper(image im, int c, int h, int w) {
    float *ptr = &im.data[0];
    int i;
    // move ptr to w
    for (i = 0; i < w; i++) {
        ptr++;
    }
    // move ptr to h
    for (i = 0; i < h; i++) {
        ptr += im.w;
    }
    // move ptr to c
    for (i = 0; i < c; i++) {
        ptr += (im.w * im.h);
    }
    return ptr;
}

// should return the pixel value at channel c, row h, and column w
// need to do bounds checking to make sure the coordinates are valid for the image
float get_pixel(image im, int c, int h, int w)
{
    // clamp bounds checking
    if (c < 0)
        c = 0;
    if (h < 0)
        h = 0;
    if (w < 0)
        w = 0;
    if (c >= im.c)
        c = im.c-1;
    if (h >= im.h)
        h = im.h-1;
    if (w >= im.w)
        w = im.w-1;

    return *get_pixel_helper(im, c, h, w);
}

// should set the pixel to the value v
// need to do bounds checking to make sure the coordinates are valid for the image
// should simply return without doing anything if you pass in invalid coordinates
void set_pixel(image im, int c, int h, int w, float v)
{
    // bounds checking
    if (c < 0 || h < 0 || w < 0)
        return;
    if (c >= im.c || h >= im.h || w >= im.w)
        return;
    
    float* ptr_to_pixel;
    ptr_to_pixel = get_pixel_helper(im, c, h, w);
    *ptr_to_pixel = v;
}

image copy_image(image im)
{
    image copy = make_image(im.c, im.h, im.w);

    int c, h, w;
    float curr_pixel;
    for (c = 0; c < im.c; c++) {
        for (h = 0; h < im.h; h++) {
            for (w = 0; w < im.w; w++) {
                curr_pixel = get_pixel(im, c, h, w);
                set_pixel(copy, c, h, w, curr_pixel);
            }
        }
    }
    return copy;
}

image rgb_to_grayscale(image im)
{
    assert(im.c == 3);
    image gray = make_image(1, im.h, im.w);
    
    int h, w;
    float c1_pixel, c2_pixel, c3_pixel, gray_pixel;
    for (h = 0; h < im.h; h++) {
        for (w = 0; w < im.w; w++) {
            // red
            c1_pixel = get_pixel(im, 0, h, w);
            // green
            c2_pixel = get_pixel(im, 1, h, w);
            // blue
            c3_pixel = get_pixel(im, 2, h, w);
            gray_pixel = 0.299*c1_pixel + 0.587*c2_pixel + 0.114*c3_pixel;
            set_pixel(gray, 0, h, w, gray_pixel);
        }
    }
    return gray;
}

// should add v to every pixel in channel c in the image
void shift_image(image im, int c, float v)
{
    int h, w;
    float curr_pixel;
    for (h = 0; h < im.h; h++) {
        for (w = 0; w < im.w; w++) {
            curr_pixel = get_pixel(im, c, h, w);
            set_pixel(im, c, h, w, curr_pixel+v);
        }
    }
}

void clamp_image(image im)
{
    int c, h, w;
    float curr_pixel;
    for (c = 0; c < im.c; c++) {
        for (h = 0; h < im.h; h++) {
            for (w = 0; w < im.w; w++) {
                curr_pixel = get_pixel(im, c, h, w);
                if (curr_pixel < 0)
                    set_pixel(im, c, h, w, 0);
                if (curr_pixel > 1)
                    set_pixel(im, c, h, w, 1);
            }
        }
    }
}

// These might be handy
float three_way_max(float a, float b, float c)
{
    return (a > b) ? ( (a > c) ? a : c) : ( (b > c) ? b : c) ;
}

float three_way_min(float a, float b, float c)
{
    return (a < b) ? ( (a < c) ? a : c) : ( (b < c) ? b : c) ;
}

// simply replace the R channel with H, the G channel with S, etc
void rgb_to_hsv(image im)
{
    int h, w;
    float hue, saturation, value, min, r, g, b;
    for (h = 0; h < im.h; h++) {
        for (w = 0; w < im.w; w++) {
            r = get_pixel(im, 0, h, w);
            g = get_pixel(im, 1, h, w);
            b = get_pixel(im, 2, h, w);

            // value is just the largest of the 3 RGB components
            value = three_way_max(r, g, b);
            // calculate saturation
            min = three_way_min(r, g, b);
            saturation = (value == 0) ? 0 : (value-min) / value;
            // calculate hue
            if ((value-min) == 0) {
                hue = 0;
            } else if (value == r) {
                hue = (g - b) / (value-min);
            } else if (value == g) {
                hue = (b - r) / (value-min) + 2.0;
            } else {
                hue = (r - g) / (value-min) + 4.0;
            }
            if (hue < 0) {
                hue = (hue / 6.0) + 1.0;
            } else if (hue > 0) {
                hue /= 6.0;
            }

            set_pixel(im, 0, h, w, hue);
            set_pixel(im, 1, h, w, saturation);
            set_pixel(im, 2, h, w, value);
        }
    }
}

void hsv_to_rgb(image im)
{
    int h, w;
    float chroma, hue, saturation, value, hue_prime, x, r, g, b, m;
    for (h = 0; h < im.h; h++) {
        for (w = 0; w < im.w; w++) {
            hue = get_pixel(im, 0, h, w);
            saturation = get_pixel(im, 1, h, w);
            value = get_pixel(im, 2, h, w);

            chroma = value * saturation;
            hue_prime = hue * 6.0;
            x = chroma * (1.0 - fabs(fmod(hue_prime, 2.0) - 1.0));
            if (0 <= hue_prime && hue_prime < 1.0) {
                r = chroma;
                g = x;
                b = 0;
            } else if (1.0 <= hue_prime && hue_prime < 2.0) {
                r = x;
                g = chroma;
                b = 0;
            } else if (2.0 <= hue_prime && hue_prime < 3.0) {
                r = 0;
                g = chroma;
                b = x;
            } else if (3.0 <= hue_prime && hue_prime < 4.0) {
                r = 0;
                g = x;
                b = chroma;
            } else if (4.0 <= hue_prime && hue_prime < 5.0) {
                r = x;
                g = 0;
                b = chroma;
            } else {
                r = chroma;
                g = 0;
                b = x;
            }

            m = value - chroma;
            r += m;
            g += m;
            b += m;
            set_pixel(im, 0, h, w, r);
            set_pixel(im, 1, h, w, g);
            set_pixel(im, 2, h, w, b);
        }
    }
}
