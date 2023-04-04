#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
#include <assert.h>
#include "image.h"
#include "matrix.h"

// Draws a line on an image with color corresponding to the direction of line
// image im: image to draw line on
// float x, y: starting point of line
// float dx, dy: vector corresponding to line angle and magnitude
void draw_line(image im, float y, float x, float dy, float dx)
{
    assert(im.c == 3);
    float angle = 6*(atan2(dy, dx) / TWOPI + .5);
    int index = floor(angle);
    float f = angle - index;
    float r, g, b;
    if(index == 0){
        r = 1; g = f; b = 0;
    } else if(index == 1){
        r = 1-f; g = 1; b = 0;
    } else if(index == 2){
        r = 0; g = 1; b = f;
    } else if(index == 3){
        r = 0; g = 1-f; b = 1;
    } else if(index == 4){
        r = f; g = 0; b = 1;
    } else {
        r = 1; g = 0; b = 1-f;
    }
    float i;
    float d = sqrt(dx*dx + dy*dy);
    for(i = 0; i < d; i += 1){
        int xi = x + dx*i/d;
        int yi = y + dy*i/d;
        set_pixel(im, 0, yi, xi, r);
        set_pixel(im, 1, yi, xi, g);
        set_pixel(im, 2, yi, xi, b);
    }
}

// dont ask me why LOL
void set_x(int c, int h, int w, float *x, image integ) {
    if (w - 1 < 0) *x = 0;
    else *x = get_pixel(integ, c, h, w-1);
}
void set_y(int c, int h, int w, float *y, image integ) {
    if (h - 1 < 0) *y = 0;
    else *y = get_pixel(integ, c, h-1, w);
}
void set_z(int c, int h, int w, float *z, image integ) {
    if ((h - 1 < 0) || (w - 1 < 0)) *z = 0;
    else *z = get_pixel(integ, c, h-1, w-1);
}

// Make an integral image or summed area table from an image
// image im: image to process
// returns: image I such that I[x,y] = sum{i<=x, j<=y}(im[i,j])
image make_integral_image(image im)
{
    image integ = copy_image(im);
    // TODO: fill in the integral image
    // x = s(w-1,h), y = s(w,h-1), z = s(w-1,h-1)
    float x, y, z;
    int h, w, c;
    float value;
    for (c = 0; c < im.c; c++) {
        for (h = 0; h < im.h; h++) {
            for (w = 0; w < im.w; w++) {
                // set s(w-1,h-1)
                set_z(c, h, w, &z, integ);
                // set s(w,h-1)
                set_y(c, h, w, &y, integ);
                // set s(w-1,h)
                set_x(c, h, w, &x, integ);
                value = get_pixel(im, c, h, w) + x + y - z;
                set_pixel(integ, c, h, w, value);
            }
        }
    }

    return integ;
}

// w & h = current position in loop over image
float box_filter_helper(image integ, int w, int h, int c, int s) {
    float q = 0.0;
    int mid = s / 2;
    float bottom_right, top_right, bottom_left, top_left;

    if (h + mid > integ.h - 1 || w + mid > integ.w - 1)
        bottom_right = 0;
    else 
        bottom_right = get_pixel(integ, c, h + mid, w + mid);

    if (h - mid - 1 < 0 || w + mid > integ.w - 1) 
        top_right = 0;
    else
        top_right = get_pixel(integ, c, h - mid - 1, w + mid);

    if (w - mid - 1 < 0 || h - mid - 1 < 0)
        top_left = 0;
    else 
        top_left = get_pixel(integ, c, h - mid - 1, w - mid - 1);

    if (h + mid > integ.h - 1 || w - mid - 1 < 0)
        bottom_left = 0;
    else 
        bottom_left = get_pixel(integ, c, h + mid, w - mid - 1);

    q = bottom_right - top_right - bottom_left + top_left;
    return q;
}

// Apply a box filter to an image using an integral image for speed
// image im: image to smooth
// int s: window size for box filter
// returns: smoothed image
image box_filter_image(image im, int s)
{
    int i,j,k;
    image integ = make_integral_image(im);
    image S = make_image(im.c, im.h, im.w);
    float q = 0.0;

    // TODO: fill in S using the integral image.
    for (k = 0; k < im.c; k++) {
        for (j = 0; j < im.h; j++) {
            for (i = 0; i < im.w; i++) {
                q = box_filter_helper(integ, i, j, k, s);
                set_pixel(S, k, j, i, q / (s*s));
            }
        }
    }

    free_image(integ); 
    return S;
}

// Calculate the time-structure matrix of an image pair.
// image im: the input image.
// image prev: the previous image in sequence.
// int s: window size for smoothing.
// returns: structure matrix. 1st channel is Ix^2, 2nd channel is Iy^2,
//          3rd channel is IxIy, 4th channel is IxIt, 5th channel is IyIt.
image time_structure_matrix(image im, image prev, int s)
{
    int converted = 0;
    if(im.c == 3){
        converted = 1;
        im = rgb_to_grayscale(im);
        prev = rgb_to_grayscale(prev);
    }

    // TODO: calculate gradients, structure components, and smooth them
    // Spatial gradients can be calculated as normal
    image S;

    image temp = make_image(5, im.h, im.w);
    // Calculate image derivatives Ix and Iy.
    // correspond to convolve im with gx & gy
    image gx_filter = make_gx_filter();
    image gy_filter = make_gy_filter();
    image gx = convolve_image(im, gx_filter, 0);
    image gy = convolve_image(im, gy_filter, 0);
    free_image(gx_filter);
    free_image(gy_filter);

    // Calculate measures IxIx, IyIy, and IxIy.
    // The time gradient can be calculated as 
    // the difference between the previous image and the next image in a sequence.
    // 4th channel is IxIt, 5th channel is IyIt.
    int h, w;
    for (h = 0; h < im.h; h++) {
        for (w = 0; w < im.w; w++) {
            float Ix = get_pixel(gx, 0, h, w);
            float Iy = get_pixel(gy, 0, h, w);
            float prev_current_pix = get_pixel(prev, 0, h, w);
            float current_pix = get_pixel(im, 0, h, w);
            set_pixel(temp, 0, h, w, powf(Ix, 2));
            set_pixel(temp, 1, h, w, powf(Iy, 2));
            set_pixel(temp, 2, h, w, Ix * Iy);
            set_pixel(temp, 3, h, w, (current_pix * Ix) - (prev_current_pix * Ix));
            set_pixel(temp, 4, h, w, (current_pix * Iy) - (prev_current_pix * Iy));
        }
    }

    // Calculate structure matrix components as weighted sum of nearby measures.
    // using cooler box filter blur B)
    S = box_filter_image(temp, s);
    free_image(temp);
    free_image(gx);
    free_image(gy);
    if(converted){
        free_image(im); free_image(prev);
    }
    return S;
}

// Calculate the velocity given a structure image
// image S: time-structure image
// int stride: 
image velocity_image(image S, int stride)
{
    image v = make_image(3, S.h/stride, S.w/stride);
    int i, j;
    matrix M = make_matrix(2,2);
    for(j = (stride-1)/2; j < S.h; j += stride){
        for(i = (stride-1)/2; i < S.w; i += stride){
            float Ixx = S.data[i + S.w*j + 0*S.w*S.h];
            float Iyy = S.data[i + S.w*j + 1*S.w*S.h];
            float Ixy = S.data[i + S.w*j + 2*S.w*S.h];
            float Ixt = S.data[i + S.w*j + 3*S.w*S.h];
            float Iyt = S.data[i + S.w*j + 4*S.w*S.h];

            // TODO: calculate vx and vy using the flow equation
            float det = Ixx*Iyy - powf(Ixy, 2);
            float vx, vy;
            // ensure that the matrix has an inverse
            if (det == 0) {
                vx = 0;
                vy = 0;
            } else {
                matrix inverse, M2, M3;
                M.data[0][0] = Ixx;
                M.data[0][1] = Ixy;
                M.data[1][0] = Ixy;
                M.data[1][1] = Iyy;
                inverse = matrix_invert(M);

                M2 = make_matrix(2, 1);
                M2.data[0][0] = -1 * Ixt;
                M2.data[1][0] = -1 * Iyt;
                M3 = matrix_mult_matrix(inverse, M2);

                vx = M3.data[0][0];
                vy = M3.data[1][0];

                free_matrix(inverse);
                free_matrix(M2);
                free_matrix(M3);
            }


            set_pixel(v, 0, j/stride, i/stride, vx);
            set_pixel(v, 1, j/stride, i/stride, vy);
        }
    }
    free_matrix(M);
    return v;
}

// Draw lines on an image given the velocity
// image im: image to draw on
// image v: velocity of each pixel
// float scale: scalar to multiply velocity by for drawing
void draw_flow(image im, image v, float scale)
{
    int stride = im.w / v.w;
    int i,j;
    for (j = (stride-1)/2; j < im.h; j += stride) {
        for (i = (stride-1)/2; i < im.w; i += stride) {
            float dx = scale*get_pixel(v, 0, j/stride, i/stride);
            float dy = scale*get_pixel(v, 1, j/stride, i/stride);
            if(fabs(dx) > im.w) dx = 0;
            if(fabs(dy) > im.h) dy = 0;
            draw_line(im, j, i, dy, dx);
        }
    }
}


// Constrain the absolute value of each image pixel
// image im: image to constrain
// float v: each pixel will be in range [-v, v]
void constrain_image(image im, float v)
{
    int i;
    for(i = 0; i < im.w*im.h*im.c; ++i){
        if (im.data[i] < -v) im.data[i] = -v;
        if (im.data[i] >  v) im.data[i] =  v;
    }
}

// Calculate the optical flow between two images
// image im: current image
// image prev: previous image
// int smooth: amount to smooth structure matrix by
// int stride: downsampling for velocity matrix
// returns: velocity matrix
image optical_flow_images(image im, image prev, int smooth, int stride)
{
    image S = time_structure_matrix(im, prev, smooth);   
    image v = velocity_image(S, stride);
    constrain_image(v, 6);
    image vs = smooth_image(v, 2);
    free_image(v);
    free_image(S);
    return vs;
}

// Run optical flow demo on webcam
// int smooth: amount to smooth structure matrix by
// int stride: downsampling for velocity matrix
// int div: downsampling factor for images from webcam
void optical_flow_webcam(int smooth, int stride, int div)
{
#ifdef OPENCV
    void * cap;
    // What video stream you open
    cap = open_video_stream(0, 1, 0, 0, 0);
    printf("%ld\n", cap);
    if(!cap){
        fprintf(stderr, "couldn't open\n");
        exit(0);
    }
    image prev = get_image_from_stream(cap);
    printf("%d %d\n", prev.w, prev.h);
    image prev_c = nn_resize(prev, prev.h/div, prev.w/div);
    image im = get_image_from_stream(cap);
    image im_c = nn_resize(im, im.h/div, im.w/div);
    while(im.data){
        image copy = copy_image(im);
        image v = optical_flow_images(im_c, prev_c, smooth, stride);
        draw_flow(copy, v, smooth*div*2);
        int key = show_image(copy, "flow", 5);
        free_image(v);
        free_image(copy);
        free_image(prev);
        free_image(prev_c);
        prev = im;
        prev_c = im_c;
        if(key != -1) {
            key = key % 256;
            printf("%d\n", key);
            if (key == 27) break;
        }
        im = get_image_from_stream(cap);
        im_c = nn_resize(im, im.h/div, im.w/div);
    }
#else
    fprintf(stderr, "Must compile with OpenCV\n");
#endif
}
