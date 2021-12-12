#include "DrawScene.h"
#include "utils/LogUtil.h"
#include <iostream>

/**
 * \brief               hard to understand, what's this?
 */
void insertImageMemoryBarrier(
    VkCommandBuffer cmdbuffer, VkImage image, VkAccessFlags srcAccessMask,
    VkAccessFlags dstAccessMask, VkImageLayout oldImageLayout,
    VkImageLayout newImageLayout, VkPipelineStageFlags srcStageMask,
    VkPipelineStageFlags dstStageMask, VkImageSubresourceRange subresourceRange)
{

    VkImageMemoryBarrier imageMemoryBarrier{};
    {
        imageMemoryBarrier.sType = VK_STRUCTURE_TYPE_IMAGE_MEMORY_BARRIER;
        imageMemoryBarrier.srcQueueFamilyIndex = VK_QUEUE_FAMILY_IGNORED;
        imageMemoryBarrier.dstQueueFamilyIndex = VK_QUEUE_FAMILY_IGNORED;
    }
    imageMemoryBarrier.srcAccessMask = srcAccessMask;
    imageMemoryBarrier.dstAccessMask = dstAccessMask;
    imageMemoryBarrier.oldLayout = oldImageLayout;
    imageMemoryBarrier.newLayout = newImageLayout;
    imageMemoryBarrier.image = image;
    imageMemoryBarrier.subresourceRange = subresourceRange;

    vkCmdPipelineBarrier(cmdbuffer, srcStageMask, dstStageMask, 0, 0, nullptr,
                         0, nullptr, 1, &imageMemoryBarrier);
}
#define DEFAULT_FENCE_TIMEOUT 100000000000
void cDrawScene::flushCommandBuffer(VkCommandBuffer commandBuffer,
                                    VkQueue queue, VkCommandPool pool,
                                    bool free)
{
    if (commandBuffer == VK_NULL_HANDLE)
    {
        return;
    }

    SIM_ASSERT(vkEndCommandBuffer(commandBuffer) == VK_SUCCESS);

    VkSubmitInfo submitInfo{};
    submitInfo.sType = VK_STRUCTURE_TYPE_SUBMIT_INFO;
    submitInfo.commandBufferCount = 1;
    submitInfo.pCommandBuffers = &commandBuffer;
    // Create fence to ensure that the command buffer has finished executing
    VkFenceCreateInfo fenceInfo{};
    fenceInfo.sType = VK_STRUCTURE_TYPE_FENCE_CREATE_INFO;
    fenceInfo.flags = 0;
    VkFence fence;
    SIM_ASSERT(vkCreateFence(mDevice, &fenceInfo, nullptr, &fence) ==
               VK_SUCCESS);
    // Submit to the queue
    SIM_ASSERT(vkQueueSubmit(queue, 1, &submitInfo, fence) == VK_SUCCESS);
    // Wait for the fence to signal that command buffer has finished executing
    SIM_ASSERT(vkWaitForFences(mDevice, 1, &fence, VK_TRUE,
                               DEFAULT_FENCE_TIMEOUT) == VK_SUCCESS);
    vkDestroyFence(mDevice, fence, nullptr);
    if (free)
    {
        vkFreeCommandBuffers(mDevice, pool, 1, &commandBuffer);
    }
}

/**
 * \brief           Take a screenshot from the color buffer
 *
 *      imitating the code from Sachas Williems
 * Take a screenshort from the current swapchain image
 * Using a blit from the swapchain image, to a linear image whose memory content
 * is then saved as a ppm image the original swapchain image cannot be saved as
 * an image directly, because they are stored in an implemention depentdent
 * optimal tilling format
 *
 * This requires the swapchain images to be created with the
 * VK_IMAGE_USAGE_TRANSFER_SRC_BIT
 */
extern uint32_t findMemoryType(VkPhysicalDevice physicalDevice,
                               uint32_t typeFilter,
                               VkMemoryPropertyFlags properties);
#include "utils/FileUtil.h"
#include <fstream>
void cDrawScene::ScreenShotDraw(std::string filename)
{
    std::string name = cFileUtil::GetDir(filename);
    // std::cout << "filename = " << filename << std::endl;
    // std::cout << "dirname = " << name << std::endl;
    if (cFileUtil::ExistsDir(name) == false)
    {
        // std::cout << "try to create dir " << name << std::endl;
        cFileUtil::CreateDir(name.c_str());
    }
    // exit(0);
    bool supportBlit = true;
    VkFormatProperties formatProps{};

    vkGetPhysicalDeviceFormatProperties(mPhysicalDevice,
                                        VK_FORMAT_B8G8R8A8_SRGB, &formatProps);

    // check if the device support blitting from optimal images (the swapchain
    // images are in optimal format)
    if (!(formatProps.optimalTilingFeatures & VK_FORMAT_FEATURE_BLIT_SRC_BIT))
    {
        std::cerr << "Device does not support blitting from optimal tilling "
                     "images, using copy instead of blit!\n";
        supportBlit = false;
    }

    // check if the device support blitting to linear images
    vkGetPhysicalDeviceFormatProperties(mPhysicalDevice,
                                        VK_FORMAT_R8G8B8A8_UNORM, &formatProps);
    if (!(formatProps.linearTilingFeatures & VK_FORMAT_FEATURE_BLIT_DST_BIT))
    {
        std::cerr << "Device does not support blitting to linear tiled images, "
                     "using copy instead of blit!\n";
        supportBlit = false;
    }

    // source for the copy is the last rendered swapchain image
    int num_of_images = mSwapChainImages.size();
    int id = num_of_images - 1;
    // printf("[screenshort] num of images %d, id %d\n", num_of_images, id);
    // std::cout << "images = " << num_of_images << std::endl;
    VkImage srcImage = mSwapChainImages[id];

    // create the linear tiled destination image to copy to and to read the
    // memory from
    VkImageCreateInfo imageCreateCI{};
    imageCreateCI.sType = VK_STRUCTURE_TYPE_IMAGE_CREATE_INFO;
    imageCreateCI.imageType = VK_IMAGE_TYPE_2D;
    imageCreateCI.format = VK_FORMAT_R8G8B8A8_UNORM;
    imageCreateCI.extent.width = mSwapChainExtent.width;
    imageCreateCI.extent.height = mSwapChainExtent.height;
    imageCreateCI.extent.depth = 1;
    imageCreateCI.arrayLayers = 1;
    imageCreateCI.mipLevels = 1;
    imageCreateCI.initialLayout = VK_IMAGE_LAYOUT_UNDEFINED;
    imageCreateCI.samples = VK_SAMPLE_COUNT_1_BIT;
    imageCreateCI.tiling = VK_IMAGE_TILING_LINEAR;
    imageCreateCI.usage = VK_IMAGE_USAGE_TRANSFER_DST_BIT;

    // create the image
    VkImage dstImage;
    SIM_ASSERT(vkCreateImage(mDevice, &imageCreateCI, nullptr, &dstImage) ==
               VK_SUCCESS);

    // create the memory to back up the image
    VkMemoryRequirements memRequirements{};
    VkMemoryAllocateInfo memAllocInfo{};
    memAllocInfo.sType = VK_STRUCTURE_TYPE_MEMORY_ALLOCATE_INFO;
    VkDeviceMemory dstImageMemory{};
    vkGetImageMemoryRequirements(mDevice, dstImage, &memRequirements);
    memAllocInfo.allocationSize = memRequirements.size;
    memAllocInfo.memoryTypeIndex =
        findMemoryType(mPhysicalDevice, memRequirements.memoryTypeBits,
                       VK_MEMORY_PROPERTY_HOST_VISIBLE_BIT |
                           VK_MEMORY_PROPERTY_HOST_COHERENT_BIT);
    SIM_ASSERT(vkAllocateMemory(mDevice, &memAllocInfo, nullptr,
                                &dstImageMemory) == VK_SUCCESS);
    SIM_ASSERT(vkBindImageMemory(mDevice, dstImage, dstImageMemory, 0) ==
               VK_SUCCESS);

    // do the actual blit from the swapchain image to our host visible
    // destination image
    VkCommandBuffer copyCmd = CreateCommandBufferTool(
        VK_COMMAND_BUFFER_LEVEL_PRIMARY, mCommandPool, true);

    // transition desitnation image to transfer destination layout

    insertImageMemoryBarrier(
        copyCmd, dstImage, 0, VK_ACCESS_TRANSFER_WRITE_BIT,
        VK_IMAGE_LAYOUT_UNDEFINED, VK_IMAGE_LAYOUT_TRANSFER_DST_OPTIMAL,
        VK_PIPELINE_STAGE_TRANSFER_BIT, VK_PIPELINE_STAGE_TRANSFER_BIT,
        VkImageSubresourceRange{VK_IMAGE_ASPECT_COLOR_BIT, 0, 1, 0, 1});
    insertImageMemoryBarrier(
        copyCmd, srcImage, VK_ACCESS_MEMORY_READ_BIT,
        VK_ACCESS_TRANSFER_READ_BIT, VK_IMAGE_LAYOUT_PRESENT_SRC_KHR,
        VK_IMAGE_LAYOUT_TRANSFER_SRC_OPTIMAL, VK_PIPELINE_STAGE_TRANSFER_BIT,
        VK_PIPELINE_STAGE_TRANSFER_BIT,
        VkImageSubresourceRange{VK_IMAGE_ASPECT_COLOR_BIT, 0, 1, 0, 1});

    // If source and destination support blit, we'll blit as this also does
    // automatic format conversion (e.g. from BGR to RGB)
    if (supportBlit)
    {
        // Define the region to blit (we will blit the whold swap chain image)
        VkOffset3D blitSize;
        blitSize.x = mSwapChainExtent.width;
        blitSize.y = mSwapChainExtent.height;
        blitSize.z = 1;

        VkImageBlit imageBlitRegion{};
        imageBlitRegion.srcSubresource.aspectMask = VK_IMAGE_ASPECT_COLOR_BIT;
        imageBlitRegion.srcSubresource.layerCount = 1;
        imageBlitRegion.srcOffsets[1] = blitSize;
        imageBlitRegion.dstSubresource.aspectMask = VK_IMAGE_ASPECT_COLOR_BIT;
        imageBlitRegion.dstSubresource.layerCount = 1;
        imageBlitRegion.dstOffsets[1] = blitSize;

        // Issue the blit commdn
        vkCmdBlitImage(copyCmd, srcImage, VK_IMAGE_LAYOUT_TRANSFER_SRC_OPTIMAL,
                       dstImage, VK_IMAGE_LAYOUT_TRANSFER_DST_OPTIMAL, 1,
                       &imageBlitRegion, VK_FILTER_NEAREST);
    }
    else
    {
        // Other use image copy (requires us to manually flip components)
        VkImageCopy imageCopyRegion{};
        imageCopyRegion.srcSubresource.aspectMask = VK_IMAGE_ASPECT_COLOR_BIT;
        imageCopyRegion.srcSubresource.layerCount = 1;
        imageCopyRegion.dstSubresource.aspectMask = VK_IMAGE_ASPECT_COLOR_BIT;
        imageCopyRegion.dstSubresource.layerCount = 1;
        imageCopyRegion.extent.width = mSwapChainExtent.width;
        imageCopyRegion.extent.height = mSwapChainExtent.height;
        imageCopyRegion.extent.depth = 1;

        vkCmdCopyImage(copyCmd, srcImage, VK_IMAGE_LAYOUT_TRANSFER_SRC_OPTIMAL,
                       dstImage, VK_IMAGE_LAYOUT_TRANSFER_DST_OPTIMAL, 1,
                       &imageCopyRegion);
    }

    // transition destination iamge to general layout, which is the required
    // layout for mapping the image memory later on
    insertImageMemoryBarrier(
        copyCmd, dstImage, VK_ACCESS_TRANSFER_WRITE_BIT,
        VK_ACCESS_MEMORY_READ_BIT, VK_IMAGE_LAYOUT_TRANSFER_DST_OPTIMAL,
        VK_IMAGE_LAYOUT_GENERAL, VK_PIPELINE_STAGE_TRANSFER_BIT,
        VK_PIPELINE_STAGE_TRANSFER_BIT,
        VkImageSubresourceRange{VK_IMAGE_ASPECT_COLOR_BIT, 0, 1, 0, 1});

    // Transition back the swap chain image after the blit is done
    insertImageMemoryBarrier(
        copyCmd, srcImage, VK_ACCESS_TRANSFER_READ_BIT,
        VK_ACCESS_MEMORY_READ_BIT, VK_IMAGE_LAYOUT_TRANSFER_SRC_OPTIMAL,
        VK_IMAGE_LAYOUT_PRESENT_SRC_KHR, VK_PIPELINE_STAGE_TRANSFER_BIT,
        VK_PIPELINE_STAGE_TRANSFER_BIT,
        VkImageSubresourceRange{VK_IMAGE_ASPECT_COLOR_BIT, 0, 1, 0, 1});

    // flush command buffer
    flushCommandBuffer(copyCmd, mGraphicsQueue, this->mCommandPool);
    int width = this->mSwapChainExtent.width, height = mSwapChainExtent.height;
    // other codes...
    {

        // Get layout of the image (including row pitch)
        VkImageSubresource subResource{VK_IMAGE_ASPECT_COLOR_BIT, 0, 0};
        VkSubresourceLayout subResourceLayout;
        vkGetImageSubresourceLayout(mDevice, dstImage, &subResource,
                                    &subResourceLayout);

        // Map image memory so we can start copying from it
        const char *data;
        vkMapMemory(mDevice, dstImageMemory, 0, VK_WHOLE_SIZE, 0,
                    (void **)&data);
        data += subResourceLayout.offset;

        std::ofstream file(filename, std::ios::out | std::ios::binary);

        // ppm header
        file << "P6\n" << width << "\n" << height << "\n" << 255 << "\n";

        // If source is BGR (destination is always RGB) and we can't use blit
        // (which does automatic conversion), we'll have to manually swizzle
        // color components
        bool colorSwizzle = false;
        // Check if source is BGR
        // Note: Not complete, only contains most common and basic BGR surface
        // formats for demonstration purposes
        if (!supportBlit)
        {
            std::vector<VkFormat> formatsBGR = {VK_FORMAT_B8G8R8A8_SRGB,
                                                VK_FORMAT_B8G8R8A8_UNORM,
                                                VK_FORMAT_B8G8R8A8_SNORM};
            colorSwizzle =
                (std::find(formatsBGR.begin(), formatsBGR.end(),
                           VK_FORMAT_B8G8R8A8_SRGB) != formatsBGR.end());
        }

        // ppm binary pixel data
        for (uint32_t y = 0; y < height; y++)
        {
            unsigned int *row = (unsigned int *)data;
            for (uint32_t x = 0; x < width; x++)
            {
                if (colorSwizzle)
                {
                    file.write((char *)row + 2, 1);
                    file.write((char *)row + 1, 1);
                    file.write((char *)row, 1);
                }
                else
                {
                    file.write((char *)row, 3);
                }
                row++;
            }
            data += subResourceLayout.rowPitch;
        }
        file.close();

        // Clean up resources
        vkUnmapMemory(mDevice, dstImageMemory);
        vkFreeMemory(mDevice, dstImageMemory, nullptr);
        vkDestroyImage(mDevice, dstImage, nullptr);

        // screenshotSaved = true;
    }
}

/**
 * \brief           a helper function to create the command buffer
 * \param
 */

VkCommandBuffer cDrawScene::CreateCommandBufferTool(VkCommandBufferLevel level,
                                                    VkCommandPool pool,
                                                    bool begin)
{
    VkCommandBufferAllocateInfo cmdBufferAllocateInfo{};
    cmdBufferAllocateInfo.sType =
        VK_STRUCTURE_TYPE_COMMAND_BUFFER_ALLOCATE_INFO;
    cmdBufferAllocateInfo.commandPool = pool;
    cmdBufferAllocateInfo.level = level;
    cmdBufferAllocateInfo.commandBufferCount = 1;
    VkCommandBuffer cmdBuffer;
    SIM_ASSERT(VK_SUCCESS == vkAllocateCommandBuffers(
                                 mDevice, &cmdBufferAllocateInfo, &cmdBuffer));

    if (begin)
    {
        VkCommandBufferBeginInfo cmdBufInfo{};
        cmdBufInfo.sType = VK_STRUCTURE_TYPE_COMMAND_BUFFER_BEGIN_INFO;
        SIM_ASSERT(VK_SUCCESS == vkBeginCommandBuffer(cmdBuffer, &cmdBufInfo));
    }
    return cmdBuffer;
}
// /**
//  * \brief           Take a screenshot from the depth buffer
// */

// extern VkFormat findDepthFormat(VkPhysicalDevice phy_device);
// int convertP6toP3(const char *fileName)
// {
//     FILE *src, *dest;
//     std::string outputFilename;
//     char magicNumber[3];
//     int height, width, depth;
//     unsigned char red, green, blue;
//     int i, j, widthCounter = 1;

//     // if (checkFileExists(fileName) == FALSE)
//     // {
//     //     printf("- Given file does not exists!\n");
//     //     return ERROR;
//     // }

//     // else
//     src = fopen(fileName, "rb");

//     // create output filename #MUST FREE ALLOCATED MEMORY#
//     outputFilename = "new.ppm";
//     // REMOVE + AFTER TESTING
//     dest = fopen(outputFilename.c_str(), "w+");

//     // check that the input file is actually in P6 format
//     fscanf(src, "%s", magicNumber);
//     fscanf(src, "\n%d %d\n%d\n", &width, &height, &depth);

//     fprintf(dest, "P3\n");
//     fprintf(dest, "#P3 converted from P6\n");
//     fprintf(dest, "%d %d\n%d\n", width, height, depth);
//     ;
//     for (i = 0; i < width * height; i++)
//     {

//         for (j = 0; j < 3; j++)
//         {
//             fread(&red, 1, 1, src);
//             fread(&green, 1, 1, src);
//             fread(&blue, 1, 1, src);
//         }

//         for (j = 0; j < 3; j++)
//             fprintf(dest, "%d %d %d ", red, green, blue);

//         if (widthCounter == width)
//         {
//             fprintf(dest, "\n");
//             widthCounter = 1;
//         }

//         else
//             widthCounter++;
//     }

//     fclose(src);
//     fclose(dest);
//     return TRUE;
// }

// void cDrawScene::ScreenShotDepth(std::string filename)
// {
//     int width = mSwapChainExtent.width, height = mSwapChainExtent.height;
//     printf("[debug] width %d height %d\n", width, height);
//     VkDeviceSize size = width * height * 4;

//     // create the destination depth buffer and its memory
//     VkBuffer dstBuffer;
//     VkDeviceMemory dstMemory;
//     CreateBuffer(size, VK_BUFFER_USAGE_TRANSFER_DST_BIT,
//                  VK_MEMORY_PROPERTY_HOST_VISIBLE_BIT |
//                  VK_MEMORY_PROPERTY_HOST_COHERENT_BIT, dstBuffer, dstMemory);
//     VkCommandBuffer copyCmd =
//     CreateCommandBufferTool(VK_COMMAND_BUFFER_LEVEL_PRIMARY, mCommandPool,
//     true);

//     VkBufferImageCopy region = {};
//     region.bufferOffset = 0;
//     region.bufferImageHeight = 0;
//     region.bufferRowLength = 0;
//     region.imageSubresource.aspectMask = VK_IMAGE_ASPECT_DEPTH_BIT;
//     region.imageSubresource.mipLevel = 0;
//     region.imageSubresource.baseArrayLayer = 0;
//     region.imageSubresource.layerCount = 1;
//     region.imageOffset = VkOffset3D{0, 0, 0};
//     region.imageExtent = VkExtent3D{mSwapChainExtent.width,
//     mSwapChainExtent.height, 1};

//     // insert barrier
//     {
//         /*

// VkCommandBuffer cmdbuffer, VkImage image, VkAccessFlags srcAccessMask,
// VkAccessFlags dstAccessMask, VkImageLayout oldImageLayout,
// VkImageLayout newImageLayout, VkPipelineStageFlags srcStageMask,
// VkPipelineStageFlags dstStageMask, VkImageSubresourceRange subresourceRange
//         */
//         insertImageMemoryBarrier(
//             copyCmd, mDepthImage,
//         );
//     }
//     vkCmdCopyImageToBuffer(
//         copyCmd,
//         mDepthImage, VK_IMAGE_LAYOUT_GENERAL,
//         dstBuffer,
//         1,
//         &region);
//     flushCommandBuffer(copyCmd, mGraphicsQueue, this->mCommandPool);

//     // map image memory
//     // std::cout << "map image done\n";
//     void *data;
//     vkMapMemory(mDevice, dstMemory, 0, size, 0, &data);

//     std::ofstream file(filename, std::ios::out | std::ios::binary);

//     // ppm header
//     file << "P6\n"
//          << width << "\n"
//          << height << "\n"
//          << 255 << "\n";

//     auto size_v = width * height;

//     // auto map = [](float f) -> uint8_t {
//     //     return (uint8_t)(f * 255.0f);
//     // };

//     {
//         tMatrixXf eigen_mat = tMatrixXf::Zero(height, width);
//         float *f_data = (float *)(data);
//         // std::cout << "sizeof(float) bytes = " << sizeof(float) <<
//         std::endl;

//         // and the value should be [0, 1]
//         for (uint32_t y = 0; y < height; y++)
//         {
//             for (uint32_t x = 0; x < width; x++)
//             {
//                 eigen_mat(y, x) = f_data[y * width + x];
//             }
//         }
//         std::ofstream output_eigen_mat("eigenmat.txt");
//         output_eigen_mat << eigen_mat << std::endl;
//         output_eigen_mat.close();
//         std::cout << "output eigen mat done\n";
//     }

//     unsigned int *row = (unsigned int *)data;
//     for (uint32_t y = 0; y < height; y++)
//     {
//         for (uint32_t x = 0; x < width; x++)
//         {
//             file.write((char *)row, 1);
//             file.write((char *)row, 1);
//             file.write((char *)row, 1);
//             // eigen_mat(y, x) = static_cast<int>(static_cast<char>((((char
//             *)row)[0])));
//             // (unsigned int)(row[0])
//             row++;
//         }
//         // data += subResourceLayout.rowPitch;
//     }
//     // std::cout << "eigen mat = \n"
//     //           << eigen_mat << std::endl;
//     file.close();

//     convertP6toP3(filename.c_str());
// }