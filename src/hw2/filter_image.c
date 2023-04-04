#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
#include <assert.h>
#include "image.h"
#define TWOPI 6.2831853

void l1_normalize(image im)
{
    float total_value = 0.0;
    int c, h, w;
    for (c = 0; c < im.c; c++) {
        for (h = 0; h < im.h; h++) {
            for (w = 0; w < im.w; w++) {
                total_value += get_pixel(im, c, h, w);
            }
        }
        for (h = 0; h < im.h; h++) {
            for (w = 0; w < im.w; w++) {
                set_pixel(im, c, h, w, get_pixel(im, c, h, w) / total_value);
            }
        }
    }
    
}

image make_box_filter(int w)
{
    image box_filter = make_image(1, w, w);
    int i, j;
    for (i = 0; i < w; i++) {
        for (j = 0; j < w; j++) {
            set_pixel(box_filter, 0, i, j, 1.0 / (float)(w*w));
        }
    }
    return box_filter;
}

// this is actually a cross-correlation
image convolve_image(image im, image filter, int preserve)
{
    // filter better have either the same number of channels as im or have 1 channel.
    assert(filter.c == im.c || filter.c == 1);

    image output;
    int h, w, c, fh, fw;
    int left_bound, right_bound, top_bound, bottom_bound, h2, w2, mid;
    float q;

    if (preserve) {
        output = make_image(im.c, im.h, im.w);
        // apply the filter to every channel & preserve
        if (filter.c == 1 && im.c != 1) {
            // loop through im
            for (c = 0; c < im.c; c++) {
                for (h = 0; h < im.h; h++) {
                    for (w = 0; w < im.w; w++) {
                        mid = filter.w / 2;
                        left_bound = w - mid;
                        right_bound = w + mid;
                        top_bound = h - mid;
                        bottom_bound = h + mid;
                        fh = 0;
                        fw = 0;
                        // loop through filter relative to curent pixel in im
                        for (h2 = top_bound, fh = 0; h2 <= bottom_bound; h2++, fh++) {
                            for (w2 = left_bound, fw = 0; w2 <= right_bound; w2++, fw++) {
                                q += get_pixel(im, c, h2, w2) * get_pixel(filter, 0, fh, fw);
                            }
                        }
                        set_pixel(output, c, h, w, q);
                        q = 0.0;
                    }
                }
            }
        } else { // else normal convolution with preserve
            // loop through im
            for (c = 0; c < im.c; c++) {
                for (h = 0; h < im.h; h++) {
                    for (w = 0; w < im.w; w++) {
                        mid = filter.w / 2;
                        left_bound = w - mid;
                        right_bound = w + mid;
                        top_bound = h - mid;
                        bottom_bound = h + mid;
                        fh = 0;
                        fw = 0;
                        // loop through filter relative to curent pixel in im
                        for (h2 = top_bound, fh = 0; h2 <= bottom_bound; h2++, fh++) {
                            for (w2 = left_bound, fw = 0; w2 <= right_bound; w2++, fw++) {
                                q += get_pixel(im, c, h2, w2) * get_pixel(filter, c, fh, fw);
                            }
                        }
                        set_pixel(output, c, h, w, q);
                        q = 0.0;
                    }
                }
            }
        }
    // no preserve - using temp image to sum them up after
    } else {
        // image temp = make_image(im.c, im.h, im.w);
        output = make_image(1, im.h, im.w);
        if (filter.c == 1 && im.c != 1) {
            // loop through im
            for (h = 0; h < im.h; h++) {
                for (w = 0; w < im.w; w++) {
                    mid = filter.w / 2;
                    left_bound = w - mid;
                    right_bound = w + mid;
                    top_bound = h - mid;
                    bottom_bound = h + mid;
                    fh = 0;
                    fw = 0;
                    for (c = 0; c < im.c; c++) {
                        // loop through filter relative to curent pixel in im
                        for (h2 = top_bound, fh = 0; h2 <= bottom_bound; h2++, fh++) {
                            for (w2 = left_bound, fw = 0; w2 <= right_bound; w2++, fw++) {
                                q += get_pixel(im, c, h2, w2) * get_pixel(filter, 0, fh, fw);
                            }
                        }
                    }
                    set_pixel(output, 0, h, w, q);
                    q = 0.0;
                }
            }
        }
        else {
            // else normal convolution without preserve
            // loop through im
            for (h = 0; h < im.h; h++) {
                for (w = 0; w < im.w; w++) {
                    mid = filter.w / 2;
                    left_bound = w - mid;
                    right_bound = w + mid;
                    top_bound = h - mid;
                    bottom_bound = h + mid;
                    fh = 0;
                    fw = 0;
                    for (c = 0; c < im.c; c++) {
                        // loop through filter relative to curent pixel in im
                        for (h2 = top_bound, fh = 0; h2 <= bottom_bound; h2++, fh++) {
                            for (w2 = left_bound, fw = 0; w2 <= right_bound; w2++, fw++) {
                                q += get_pixel(im, c, h2, w2) * get_pixel(filter, c, fh, fw);
                            }
                        }
                    }
                    set_pixel(output, 0, h, w, q);
                    q = 0.0;
                }
            }
        }
    }

    return output;
}

image make_highpass_filter()
{
    image highpass = make_image(1, 3, 3);
    // corners          c, h, w, v
    set_pixel(highpass, 0, 0, 0, 0);
    set_pixel(highpass, 0, 0, 2, 0);
    set_pixel(highpass, 0, 2, 0, 0);
    set_pixel(highpass, 0, 2, 2, 0);
    // outer mids
    set_pixel(highpass, 0, 0, 1, -1);
    set_pixel(highpass, 0, 2, 1, -1);
    set_pixel(highpass, 0, 1, 0, -1);
    set_pixel(highpass, 0, 1, 2, -1);
    // middle
    set_pixel(highpass, 0, 1, 1, 4);

    return highpass;
}

image make_sharpen_filter()
{
    image sharpen = make_image(1, 3, 3);
    // corners         c, h, w, v
    set_pixel(sharpen, 0, 0, 0, 0);
    set_pixel(sharpen, 0, 0, 2, 0);
    set_pixel(sharpen, 0, 2, 0, 0);
    set_pixel(sharpen, 0, 2, 2, 0);
    // outer mids
    set_pixel(sharpen, 0, 0, 1, -1);
    set_pixel(sharpen, 0, 2, 1, -1);
    set_pixel(sharpen, 0, 1, 0, -1);
    set_pixel(sharpen, 0, 1, 2, -1);
    // middle
    set_pixel(sharpen, 0, 1, 1, 5);

    return sharpen;
}

image make_emboss_filter()
{
    image emboss = make_image(1, 3, 3);

    set_pixel(emboss, 0, 0, 0, -2);
    set_pixel(emboss, 0, 0, 2, 0);
    set_pixel(emboss, 0, 2, 0, 0);
    set_pixel(emboss, 0, 2, 2, 2);
    // outer mids
    set_pixel(emboss, 0, 0, 1, -1);
    set_pixel(emboss, 0, 2, 1, 1);
    set_pixel(emboss, 0, 1, 0, -1);
    set_pixel(emboss, 0, 1, 2, 1);
    // middle
    set_pixel(emboss, 0, 1, 1, 1);


    return emboss;
}

// Question 2.2.1: Which of these filters should we use preserve when we run our convolution and which ones should we not? Why?
// Answer: We should preserve for all when using these filters for convolutions because if we don't then our 
    // filter might impact our images' colors in a way we don't necessarily want. If I had to guess though, I would say that we should DEFINITELY preserve 
    // when running sharpen & emboss, because these filters enahance colors. Maybe we dont necessarily need to preserve when using high pass because 
    // the this filter is typically used for edge detection and not for making colors look nice.

// Question 2.2.2: Do we have to do any post-processing for the above filters? Which ones and why?
// Answer: Yes, sharpen, high pass, and emboss all need post-processing applied. To be specific after running convolutions with the image we need to 
    // clamp the image because some pixel values might have gone out of bounds (outside of the range 0 - 1.0).

// take a standard deviation value and return a filter that smooths using a gaussian with that sigma
image make_gaussian_filter(float sigma)
{
    image guassian;
    int dimensions, h, w;
    float q;
    dimensions = ((int)sigma * 6 % 2 == 1) ? sigma * 6 : sigma * 6 + 1;
    guassian = make_image(1, dimensions, dimensions);

    for (h = 0; h < dimensions; h++) {
        for (w = 0; w < dimensions; w++) {
            float x = h - (dimensions - 1) / 2.0;
            float y = w - (dimensions - 1) / 2.0;
            q = 1.0 / (TWOPI * powf(sigma, 2)) * expf(-((powf(x, 2) + powf(y, 2)) / ((2 * powf(sigma, 2)))));
            set_pixel(guassian, 0, h, w, q);
        }
    }

    l1_normalize(guassian);
    return guassian;
}

image add_image(image a, image b)
{
    assert(a.h == b.h && a.w == b.w && a.c == b.c);

    int h, w, c;
    float q = 0;
    image sum_im = make_image(a.c, a.h, a.w);

    for (c = 0; c < a.c; c++) {
        for (h = 0; h < a.h; h++) {
            for (w = 0; w < a.w; w++) {
                q = get_pixel(a, c, h, w) + get_pixel(b, c, h, w);
                set_pixel(sum_im, c, h, w, q);
            }
        }
    }
    return sum_im;
}

image sub_image(image a, image b)
{
    assert(a.h == b.h && a.w == b.w && a.c == b.c);

    int h, w, c;
    float q = 0;
    image sub_im = make_image(a.c, a.h, a.w);

    for (c = 0; c < a.c; c++) {
        for (h = 0; h < a.h; h++) {
            for (w = 0; w < a.w; w++) {
                q = get_pixel(a, c, h, w) - get_pixel(b, c, h, w);
                set_pixel(sub_im, c, h, w, q);
            }
        }
    }
    return sub_im;
}

image make_gx_filter()
{
    image gx = make_image(1, 3, 3);
    // corners    c, h, w, v
    set_pixel(gx, 0, 0, 0, -1);
    set_pixel(gx, 0, 0, 2, 1);
    set_pixel(gx, 0, 2, 0, -1);
    set_pixel(gx, 0, 2, 2, 1);
    // outer mids
    set_pixel(gx, 0, 0, 1, 0);
    set_pixel(gx, 0, 2, 1, 0);
    set_pixel(gx, 0, 1, 0, -2);
    set_pixel(gx, 0, 1, 2, 2);
    // middle
    set_pixel(gx, 0, 1, 1, 0);

    return gx;
}

image make_gy_filter()
{
    image gy = make_image(1, 3, 3);
    // corners    c, h, w, v
    set_pixel(gy, 0, 0, 0, -1);
    set_pixel(gy, 0, 0, 2, -1);
    set_pixel(gy, 0, 2, 0, 1);
    set_pixel(gy, 0, 2, 2, 1);
    // outer mids
    set_pixel(gy, 0, 0, 1, -2);
    set_pixel(gy, 0, 2, 1, 2);
    set_pixel(gy, 0, 1, 0, 0);
    set_pixel(gy, 0, 1, 2, 0);
    // middle
    set_pixel(gy, 0, 1, 1, 0);

    return gy;
}

void feature_normalize(image im)
{
    float curr_pixel, normalized;
    int c, h, w;

    float min = im.data[0];
    float max = im.data[0];
    for (c = 0; c < im.c; c++) {
        for (h = 0; h < im.h; h++) {
            for (w = 0; w < im.w; w++) {
                curr_pixel = get_pixel(im, c, h, w);
                if (curr_pixel < min)
                    min = curr_pixel;
                if (curr_pixel > max)
                    max = curr_pixel;
            }
        }
    }

    // scale the image so all values lie between [0-1]
    for (c = 0; c < im.c; c++) {
        for (h = 0; h < im.h; h++) {
            for (w = 0; w < im.w; w++) {
                curr_pixel = get_pixel(im, c, h, w);
                normalized = (curr_pixel - min) / (max - min);
                set_pixel(im, c, h, w, normalized);
            }
        }
    }

}

image *sobel_image(image im)
{
    image* im2 = calloc(2, sizeof(image));

    // convolve im with gx & gy
    image gx_filter = make_gx_filter();
    image gy_filter = make_gy_filter();
    image gx = convolve_image(im, gx_filter, 0);
    image gy = convolve_image(im, gy_filter, 0);

    // get gradient magnitude
    int h, w, h_max, w_max;
    float q;
    h_max = im.h;
    w_max = im.w;
    image mag = make_image(1, im.h, im.w);
    for (h = 0; h < h_max; h++) {
        for (w = 0; w < w_max; w++) {
            q = sqrtf(powf(get_pixel(gx, 0, h, w), 2) + powf(get_pixel(gy, 0, h, w), 2));
            set_pixel(mag, 0, h, w, q);
        }
    }

    // get gradient's direction (theta)
    image theta = make_image(1, im.h, im.w);
    for (h = 0; h < h_max; h++) {
        for (w = 0; w < w_max; w++) {
            q = atan2f(get_pixel(gy, 0, h, w), get_pixel(gx, 0, h, w));
            set_pixel(theta, 0, h, w, q);
        }
    }

    free_image(gx_filter);
    free_image(gy_filter);
    free_image(gx);
    free_image(gy);
    
    im2[0] = mag;
    im2[1] = theta;

    return im2;
}

image colorize_sobel(image im)
{
    // TODO
    return make_image(1,1,1);
}
