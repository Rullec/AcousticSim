#define TINYEXR_IMPLEMENTATION
#include "tinyexr.h"

/**
 * \brief           save depth exr
 */
bool SaveEXRSingleChannel(const float *depth, int width, int height,
                          const char *outfilename)
{
    EXRHeader header;
    InitEXRHeader(&header);
    header.compression_type = TINYEXR_COMPRESSIONTYPE_PIZ;

    EXRImage image;
    InitEXRImage(&image);

    image.num_channels = 1;

    std::vector<float> images[1];
    images[0].resize(width * height);
    // images[1].resize(width * height);
    // images[2].resize(width * height);

    // Split RGBRGBRGB... into R, G and B layer
    for (int i = 0; i < width * height; i++)
    {
        images[0][i] = depth[i + 0];
        // images[1][i] = depth[3 * i + 1];
        // images[2][i] = depth[3 * i + 2];
    }

    float *image_ptr[3];
    image_ptr[0] = &(images[0].at(0)); // B
    // image_ptr[1] = &(images[1].at(0)); // G
    // image_ptr[2] = &(images[0].at(0)); // R

    image.images = (unsigned char **)image_ptr;
    image.width = width;
    image.height = height;

    header.num_channels = 1;
    header.channels =
        (EXRChannelInfo *)malloc(sizeof(EXRChannelInfo) * header.num_channels);
    // Must be (A)BGR order, since most of EXR viewers expect this channel
    // order. strncpy(header.channels[0].name, "B", 255);
    // header.channels[0].name[strlen("B")] = '\0';
    // strncpy(header.channels[1].name, "G", 255);
    // header.channels[1].name[strlen("G")] = '\0';
    strncpy(header.channels[0].name, "Z", 255);
    header.channels[0].name[strlen("Z")] = '\0';

    header.pixel_types = (int *)malloc(sizeof(int) * header.num_channels);
    header.requested_pixel_types =
        (int *)malloc(sizeof(int) * header.num_channels);
    for (int i = 0; i < header.num_channels; i++)
    {
        header.pixel_types[i] =
            TINYEXR_PIXELTYPE_FLOAT; // pixel type of input image
        // header.requested_pixel_types[i] = TINYEXR_PIXELTYPE_HALF; // pixel
        // type of output image to be stored in .EXR
        header.requested_pixel_types[i] =
            TINYEXR_PIXELTYPE_FLOAT; // pixel type of output image to be stored
                                     // in .EXR
    }

    const char *err = NULL; // or nullptr in C++11 or later.
    int ret = SaveEXRImageToFile(&image, &header, outfilename, &err);
    if (ret != TINYEXR_SUCCESS)
    {
        fprintf(stderr, "Save EXR err: %s\n", err);
        FreeEXRErrorMessage(err); // free's buffer for an error message
        return ret;
    }
    printf("Saved exr file 1. [ %s ] \n", outfilename);

    // free(depth);

    free(header.channels);
    free(header.pixel_types);
    free(header.requested_pixel_types);
    return true;
}

/**
 * \brief               Save png images
 */

#define STB_IMAGE_WRITE_IMPLEMENTATION
#include "utils/stb_image_write.h"
#include <iostream>
bool SavePNGSingleChannel(const float *depth_pixels, int width, int height,
                          const char *outfile_name)
{
    int num_of_pixels = width * height;
    uint8_t *png_pixels = new uint8_t[num_of_pixels];
    uint8_t min_val = 255, max_val = 0;
    for (int i = 0; i < num_of_pixels; i++)
    {
        // const int int_v = int(depth_pixels[i] * 255.99f);
        png_pixels[i] = uint8_t(depth_pixels[i] * 255.99f);
        if (png_pixels[i] < min_val)
            min_val = png_pixels[i];
        if (max_val < png_pixels[i])
            max_val = png_pixels[i];

        // = rgba;
    }
    // printf("max pixel %d min pixel %d\n", max_val, min_val);
    stbi_write_png(outfile_name, width, height, 1, png_pixels,
                   width * sizeof(uint8_t));
    // std::cout << "[debug] save png image to " << outfile_name << std::endl;
    delete[] png_pixels;
    return true;
}
#include "utils/LogUtil.h"
bool SavePNGSingleChannelDepth(const float *depth_pixels, int width, int height,
                          const char *outfile_name, float range_st_m,
                          float range_ed_m)
{
    SIM_ASSERT(range_st_m < range_ed_m);
    int num_of_pixels = width * height;
    uint8_t *png_pixels = new uint8_t[num_of_pixels];
    uint8_t min_val = 255, max_val = 0;
    for (int i = 0; i < num_of_pixels; i++)
    {
        float cur_depth_val = depth_pixels[i];
        if ((cur_depth_val < range_st_m) || (cur_depth_val > range_ed_m))
        {
            SIM_ERROR("cur depth val {} out of range [{}, {}]", cur_depth_val,
                      range_st_m, range_ed_m);
            exit(1);
        }
        cur_depth_val =
            (cur_depth_val - range_st_m) / (range_ed_m - range_st_m);
        png_pixels[i] = uint8_t(cur_depth_val * 255.99f);
        if (png_pixels[i] < min_val)
            min_val = png_pixels[i];
        if (max_val < png_pixels[i])
            max_val = png_pixels[i];
    }
    // printf("max pixel %d min pixel %d\n", max_val, min_val);
    stbi_write_png(outfile_name, width, height, 1, png_pixels,
                   width * sizeof(uint8_t));
    // std::cout << "[debug] save png image to " << outfile_name << std::endl;
    delete[] png_pixels;
    return true;
}