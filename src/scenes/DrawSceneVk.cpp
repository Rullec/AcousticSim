#include "DrawScene.h"
#include "scenes/SimScene.h"
#include "utils/LogUtil.h"
#include "utils/MathUtil.h"
#include "vulkan/vulkan.h"
#include <iostream>
#include <optional>
#include <set>
#define GLM_FORCE_RADIANS
#define GLM_FORCE_DEPTH_ZERO_TO_ONE
#define GLFW_INCLUDE_VULKAN
#ifdef _WIN32
#define VK_USE_PLATFORM_WIN32_KHR
#include "GLFW/glfw3.h"
#define GLFW_EXPOSE_NATIVE_WIN32
#include <GLFW/glfw3native.h>
#endif

#ifdef __linux__
#define VK_USE_PLATFORM_XCB_KHR
#include <GLFW/glfw3.h>
#define GLFW_EXPOSE_NATIVE_X11
#include <GLFW/glfw3native.h>
#endif

#ifdef __APPLE__
#include <GLFW/glfw3.h>
#endif

extern std::vector<const char *> validationLayers;
extern std::vector<const char *> deviceExtensions;
extern GLFWwindow *window;
extern bool enableValidationLayers;

/**
 * \brief           the details about the support for swapchain, on current
 * surface and physical device
 */
struct SwapChainSupportDetails
{
    VkSurfaceCapabilitiesKHR capabilities;
    std::vector<VkSurfaceFormatKHR> formats;
    std::vector<VkPresentModeKHR> presentModes;
};

/**
 * \brief           check whehter swap chain is supported and how it's supported
 */
SwapChainSupportDetails querySwapChainSupport(VkPhysicalDevice device,
                                              VkSurfaceKHR surface)
{
    SwapChainSupportDetails details;

    // 1. capabilities
    vkGetPhysicalDeviceSurfaceCapabilitiesKHR(device, surface,
                                              &details.capabilities);

    // 2. formats
    uint32_t formatCount;
    vkGetPhysicalDeviceSurfaceFormatsKHR(device, surface, &formatCount,
                                         nullptr);
    details.formats.resize(formatCount);
    vkGetPhysicalDeviceSurfaceFormatsKHR(device, surface, &formatCount,
                                         details.formats.data());

    // 3. presentModes
    uint32_t modesCount;
    vkGetPhysicalDeviceSurfacePresentModesKHR(device, surface, &modesCount,
                                              nullptr);

    details.presentModes.resize(modesCount);
    vkGetPhysicalDeviceSurfacePresentModesKHR(device, surface, &modesCount,
                                              details.presentModes.data());

    return details;
}

/**
 * \brief       Given an physical device, find whehter all queue families we
 * want is in it If an queue family is supported, it has an value otherwise, its
 * value is <optional>::novalue
 */
QueueFamilyIndices findQueueFamilies(VkPhysicalDevice device,
                                     VkSurfaceKHR surface)
{
    QueueFamilyIndices indices;
    uint32_t num_of_families = 0;

    // 1. get all queue faimilies on this physical device
    vkGetPhysicalDeviceQueueFamilyProperties(device, &num_of_families, nullptr);
    std::vector<VkQueueFamilyProperties> families(
        num_of_families); // storage the properties of each queue on this
                          // physical device. such as maximium number of queues,
                          // types of queues...

    vkGetPhysicalDeviceQueueFamilyProperties(device, &num_of_families,
                                             families.data());

    // 2. find a queue to support graphics
    for (int i = 0; i < families.size(); i++)
    {
        const auto &x = families[i];
        // std::cout << i << " queue flags = " << x.queueFlags << std::endl;
        if (x.queueFlags & VK_QUEUE_GRAPHICS_BIT)
        {
            indices.graphicsFamily = i;
        }

        VkBool32 presentSupport = false;
        vkGetPhysicalDeviceSurfaceSupportKHR(device, i, surface,
                                             &presentSupport);
        // &presentSupport);
        if (presentSupport == true)
        {
            indices.presentFamily = i;
        }

        // stop when we get all we want
        if (indices.IsComplete() == true)
            break;
    }

    return indices;
}

/**
 * \brief       select the best surface format(color channels and types)
 */
VkSurfaceFormatKHR
chooseSwapSurfaceFormat(const std::vector<VkSurfaceFormatKHR> availableFormats)
{
    for (auto &x : availableFormats)
    {
        if (x.format == VK_FORMAT_B8G8R8A8_SRGB &&
            x.colorSpace == VK_COLOR_SPACE_SRGB_NONLINEAR_KHR)
            return x;
    }
    return availableFormats[0];
}

/**
 * \brief       choose the swapchain present mode
 *
 *          the author thought the "mailbox" which can support the triple
 * buffering is the best option
 */
VkPresentModeKHR
chooseSwapPresentKHR(const std::vector<VkPresentModeKHR> &modes)
{
    for (const auto &x : modes)
    {
        if (x == VK_PRESENT_MODE_MAILBOX_KHR)
        {
            return x;
        }
    }
    return VK_PRESENT_MODE_FIFO_KHR;
}

/**
 * \brief       choose the swap extent (framebuffer size)
 *
 */
VkExtent2D chooseSwapExtent(const VkSurfaceCapabilitiesKHR &capabilities)
{
    if (capabilities.currentExtent.width != UINT32_MAX)
    {
        return capabilities.currentExtent;
    }
    else
    {
        int width, height;
        glfwGetFramebufferSize(window, &width, &height);

        VkExtent2D actualExtent = {static_cast<uint32_t>(width),
                                   static_cast<uint32_t>(height)};

        actualExtent.width = std::max(
            capabilities.minImageExtent.width,
            std::min(capabilities.maxImageExtent.width, actualExtent.width));
        actualExtent.height = std::max(
            capabilities.minImageExtent.height,
            std::min(capabilities.maxImageExtent.height, actualExtent.height));

        return actualExtent;
    }
}
extern bool checkDeviceExtensionSupport(VkPhysicalDevice device);
/**
 * \brief       Check whether an device is suitable for us
 *      Devices are not equal
 */
bool IsDeviceSuitable(VkPhysicalDevice &dev, VkSurfaceKHR surface)
{
    // VkPhysicalDeviceProperties props;
    // VkPhysicalDeviceFeatures feas;
    // vkGetPhysicalDeviceProperties(dev, &props);

    // vkGetPhysicalDeviceFeatures(dev, &feas);

    // // discrete GPU and geo shader support?
    // return props.deviceType == VK_PHYSICAL_DEVICE_TYPE_DISCRETE_GPU &&
    // feas.geometryShader; return true; is graphics queue (or any queue we
    // want) families supported on this dev?
    QueueFamilyIndices indices = findQueueFamilies(dev, surface);

    // if some extensions are not supported, it cannot be a suitable device
    bool extensionsSupported = checkDeviceExtensionSupport(dev);

    bool swapChainAdequate = false;
    if (extensionsSupported)
    {
        auto details = querySwapChainSupport(dev, surface);
        swapChainAdequate = details.formats.empty() == false &&
                            (details.presentModes.empty()) == false;
    }
    VkPhysicalDeviceFeatures supportedFeatures;
    vkGetPhysicalDeviceFeatures(dev, &supportedFeatures);
    // std::cout << "indices complete " << indices.IsComplete() << std::endl;
    // std::cout << "extensions support " << extensionsSupported << std::endl;
    // std::cout << "swapChainAdequate " << swapChainAdequate << std::endl;
    // std::cout << "supportedFeatures.samplerAnisotropy " <<
    // supportedFeatures.samplerAnisotropy << std::endl;
    return indices.IsComplete() && extensionsSupported && swapChainAdequate &&
           supportedFeatures.samplerAnisotropy;
}

/**
 * \brief       rating a physical device
 */
int RateDeviceSutability(VkPhysicalDevice &dev, VkSurfaceKHR surface)
{
    int score = 0;
    VkPhysicalDeviceProperties props;
    VkPhysicalDeviceFeatures feas;
    vkGetPhysicalDeviceProperties(dev, &props);

    vkGetPhysicalDeviceFeatures(dev, &feas);

    if (IsDeviceSuitable(dev, surface) == false)
    {
        std::cout << "device is not suitable\n";
        return 0;
    }
    else
    {
        score += 1;
    }
    if (props.deviceType == VK_PHYSICAL_DEVICE_TYPE_DISCRETE_GPU)
        score += 1e3;
    // if (feas.geometryShader == false)
    //     return 0;
    return score;
}

VkSampleCountFlagBits GetMaxUsableSampleCount(VkPhysicalDevice device)
{
    VkPhysicalDeviceProperties physicalDeviceProperties;
    vkGetPhysicalDeviceProperties(device, &physicalDeviceProperties);

    // both depth and color
    VkSampleCountFlags counts =
        physicalDeviceProperties.limits.framebufferColorSampleCounts &
        physicalDeviceProperties.limits.framebufferDepthSampleCounts;

    if (counts & VK_SAMPLE_COUNT_64_BIT)
        return VK_SAMPLE_COUNT_64_BIT;
    if (counts & VK_SAMPLE_COUNT_32_BIT)
        return VK_SAMPLE_COUNT_32_BIT;
    if (counts & VK_SAMPLE_COUNT_16_BIT)
        return VK_SAMPLE_COUNT_16_BIT;
    if (counts & VK_SAMPLE_COUNT_8_BIT)
        return VK_SAMPLE_COUNT_8_BIT;
    if (counts & VK_SAMPLE_COUNT_4_BIT)
        return VK_SAMPLE_COUNT_4_BIT;
    if (counts & VK_SAMPLE_COUNT_2_BIT)
        return VK_SAMPLE_COUNT_2_BIT;

    return VK_SAMPLE_COUNT_1_BIT;
}

/**
 * \brief           Create Vulkan Instance
 *      create the application appInfo struct and write down some value
 */
void cDrawScene::CreateInstance()
{
    // 1. check validation layer enable & supported?
    if (enableValidationLayers == true &&
        CheckValidationLayerSupport() == false)
    {
        SIM_ERROR("Validation Layers is not supported");
    }
    // 2. application info
    VkApplicationInfo
        appInfo{}; // it contains a pNext member for further extension
    appInfo.sType =
        VK_STRUCTURE_TYPE_APPLICATION_INFO; // type should be set manually
    appInfo.pApplicationName = "AcousticSim";
    appInfo.applicationVersion = VK_MAKE_VERSION(1, 0, 0);
    appInfo.pEngineName = "No Engine";
    appInfo.engineVersion = VK_MAKE_VERSION(1, 0, 0);
    appInfo.apiVersion = VK_API_VERSION_1_0;

    // Lots' of information are passed by structure instead of function
    // parameters tell the Vulkan driver whcih global extensions and vlidation
    // layers we want to use it effects the entire program.
    VkInstanceCreateInfo createInfo{};
    createInfo.sType = VK_STRUCTURE_TYPE_INSTANCE_CREATE_INFO; // create appInfo
    createInfo.pApplicationInfo = &appInfo;
    // createInfo.pApplicationInfo = nullptr;

    {
        // vulkan is platform agnostic, we need to know what extension glfw
        // init, in order to pass it to the vulkan
        uint32_t glfwExtensionCount = 0;
        const char **glfwExtensions;
        glfwExtensions = glfwGetRequiredInstanceExtensions(&glfwExtensionCount);

        uint32_t total_extension_count = glfwExtensionCount + 2;
        const char **total_extensions = new const char *[total_extension_count];
        for (int i = 0; i < glfwExtensionCount; i++)
        {
            total_extensions[i] = glfwExtensions[i];
        }
        total_extensions[total_extension_count - 1] =
            "VK_KHR_get_physical_device_properties2";
        total_extensions[total_extension_count - 2] = "VK_EXT_debug_report";
        createInfo.enabledExtensionCount = total_extension_count;
        createInfo.ppEnabledExtensionNames = total_extensions;
        // for (int i = 0; i < total_extension_count; i++)
        // {
        //     std::cout << total_extensions[i] << std::endl;
        // }
        // //
        // exit(0);
    }

    if (enableValidationLayers == true)
    {
        // how many layers do we want to enable?
        createInfo.enabledLayerCount = validationLayers.size();
        // the layer name
        createInfo.ppEnabledLayerNames = validationLayers.data();
    }
    else
    {
        createInfo.enabledLayerCount = 0;
    }

    // given the instance info, create the instance
    SIM_ASSERT(vkCreateInstance(&createInfo, nullptr, &mInstance) ==
               VK_SUCCESS);
    // SIM_INFO("CreateInstance succ");
    // CheckAvaliableExtensions();
    // exit(0);
}

/**
 * \brief           Check avaiable extensions
 */
void cDrawScene::CheckAvaliableExtensions() const
{
    uint32_t size;
    vkEnumerateInstanceExtensionProperties(nullptr, &size, nullptr);

    std::vector<VkExtensionProperties> extension(size);
    vkEnumerateInstanceExtensionProperties(nullptr, &size, extension.data());
    for (auto &x : extension)
    {
        std::cout << x.extensionName << " supported\n";
    }
    // exit(0);
}

/**
 * \brief           Check whether validation layer is supported
 */
bool cDrawScene::CheckValidationLayerSupport() const
{
    // 1. get all supported layers
    uint32_t num_of_layer;
    vkEnumerateInstanceLayerProperties(&num_of_layer, nullptr);

    std::vector<VkLayerProperties> props(num_of_layer);
    vkEnumerateInstanceLayerProperties(&num_of_layer, props.data());

    // 2. check whether the validation layer is in it
    bool supported = false;
    for (auto &support_name : props)
    {
        // std::cout << "supported = " << support_name.layerName << std::endl;
        for (auto &requrested : validationLayers)
        {

            if (strcmp(support_name.layerName, requrested) == 0)
            {
                supported = true;
                // std::cout << "validation supported";
                break;
            }
        }
    }

    // 3. return the result
    return supported;
    // return false;
}

/**
 * \brief           set the messege callback for validation layer (not impl yet)
 */
void cDrawScene::SetupDebugMessenger() {}

/**
 * \brief           After the initializetion of device, create and select a
 * physical device
 *
 *      the physical device will be stored in the `VkPhysicalDevice`. they will
 * be released along with the `VkDevice`, do not need to release manually
 */
#include <map>
void cDrawScene::PickPhysicalDevice()
{
    mPhysicalDevice = VK_NULL_HANDLE;

    // get the number of physical devices
    uint32_t deviceCount = 0;
    vkEnumeratePhysicalDevices(mInstance, &deviceCount, nullptr);
    SIM_ASSERT(deviceCount != 0);

    // get all physical devices
    std::vector<VkPhysicalDevice> devices(deviceCount);
    vkEnumeratePhysicalDevices(mInstance, &deviceCount, devices.data());

    // rate all phy devices
    std::multimap<int, VkPhysicalDevice> candidates;
    for (int i = 0; i < devices.size(); i++)
    {
        auto x = devices[i];
        int score = RateDeviceSutability(x, mSurface);
        // std::cout << "device " << i << " score " << score << std::endl;
        candidates.insert(std::make_pair(RateDeviceSutability(x, mSurface), x));
    }
    SIM_ASSERT(candidates.size() > 0);
    // fetch the highest score physical device
    SIM_ASSERT(candidates.begin()->first > 0);
    // std::cout << candidates.begin()->first << std::endl;
    mPhysicalDevice = candidates.begin()->second;

    // get the max sample count of the device
    mSampleCount = GetMaxUsableSampleCount(mPhysicalDevice);
    // std::cout << "[debug] pick up physical device done, device =  " <<
    // mPhysicalDevice << " score = " << candidates.begin()->first << std::endl;
    SIM_ASSERT(mPhysicalDevice != VK_NULL_HANDLE);
}

/**
 * \brief           Create the logical device from the given physical device
 *
 *      write the member "VkDevice mDevice"
 */

void cDrawScene::CreateLogicalDevice()
{
    // 1. fetch all supported queue families from the physical device
    QueueFamilyIndices indices = findQueueFamilies(mPhysicalDevice, mSurface);
    // std::cout << "physical score = " << RateDeviceSutability(mPhysicalDevice)
    // << std::endl; std::cout << "has value = " <<
    // indices.graphicsFamily.has_value() << std::endl; exit(0);

    // 2. set the queue family info into the physical device
    std::vector<VkDeviceQueueCreateInfo> queueCreateinfos;
    std::set<uint32_t> uniqueQueueFamilies = {indices.graphicsFamily.value(),
                                              indices.presentFamily.value()};

    float queuePriority = 1.0f;
    for (uint32_t queueFamily : uniqueQueueFamilies)
    {
        VkDeviceQueueCreateInfo queueCreateInfo{};
        queueCreateInfo.sType = VK_STRUCTURE_TYPE_DEVICE_QUEUE_CREATE_INFO;
        queueCreateInfo.queueFamilyIndex = queueFamily;
        queueCreateInfo.queueCount = 1;
        queueCreateInfo.pQueuePriorities = &queuePriority;
        queueCreateinfos.push_back(queueCreateInfo);
    }

    // 3.
    VkPhysicalDeviceFeatures deviceFeatures{};
    deviceFeatures.samplerAnisotropy = VK_TRUE;
    deviceFeatures.sampleRateShading =
        VK_TRUE; // enable sampling shading for interior texture anti-aliasing

    // 4. create logical device
    VkDeviceCreateInfo createInfo{};
    createInfo.sType = VK_STRUCTURE_TYPE_DEVICE_CREATE_INFO;
    createInfo.pQueueCreateInfos = queueCreateinfos.data();
    createInfo.queueCreateInfoCount = queueCreateinfos.size();
    createInfo.pEnabledFeatures = &deviceFeatures;

    {
        // we don't need any more extensions here now
        {
            // set extensions count in the DeviceCreateInfo
            createInfo.enabledExtensionCount =
                static_cast<uint32_t>(deviceExtensions.size());
            createInfo.ppEnabledExtensionNames = deviceExtensions.data();
        }

        if (enableValidationLayers == true)
        {
            // how many layers do we want to enable?
            createInfo.enabledLayerCount = validationLayers.size();
            // the layer name
            createInfo.ppEnabledLayerNames = validationLayers.data();
        }
        else
        {
            createInfo.enabledLayerCount = 0;
        }
    }

    if (vkCreateDevice(mPhysicalDevice, &createInfo, nullptr, &mDevice) !=
        VK_SUCCESS)
    {
        SIM_ERROR("create logic device failed");
    }

    // get the queue handles: graphics queue and presentQueue
    vkGetDeviceQueue(mDevice, indices.graphicsFamily.value(), 0,
                     &mGraphicsQueue);
    vkGetDeviceQueue(mDevice, indices.presentFamily.value(), 0, &mPresentQueue);
}

/**
 * \brief       Create window surface for display
 */

void cDrawScene::CreateSurface()
{
    if (VK_SUCCESS !=
        glfwCreateWindowSurface(mInstance, window, nullptr, &mSurface))
    {
        SIM_ERROR("glfw create surface for vulkan failed");
    }
}

/**
 * \brief       create the swap chain
 */
void cDrawScene::CreateSwapChain()
{
    // std::cout << "begin to create swap chain\n";
    // 1. fetch the details from physical device
    auto details = querySwapChainSupport(mPhysicalDevice, mSurface);
    // 2. set the format
    VkSurfaceFormatKHR format = chooseSwapSurfaceFormat(details.formats);
    // 3. set the present mode
    VkPresentModeKHR mode = chooseSwapPresentKHR(details.presentModes);
    // 4. set the extent (resolution of buffers)
    VkExtent2D extent = chooseSwapExtent(details.capabilities);

    // 5. set the size of swap chain
    uint32_t image_count = details.capabilities.minImageCount + 1;
    image_count = std::min(details.capabilities.maxImageCount, image_count);

    // 6. set up a big structure
    VkSwapchainCreateInfoKHR create_info{};
    create_info.sType = VK_STRUCTURE_TYPE_SWAPCHAIN_CREATE_INFO_KHR;
    create_info.surface = mSurface;

    create_info.minImageCount = image_count;
    create_info.imageFormat = format.format;
    create_info.imageColorSpace = format.colorSpace;
    create_info.imageExtent = extent;
    create_info.imageArrayLayers = 1; // always 1

    // we also use the swapchan image as the source images when screenshoting
    create_info.imageUsage =
        VK_IMAGE_USAGE_COLOR_ATTACHMENT_BIT |
        VK_IMAGE_USAGE_TRANSFER_SRC_BIT; // we use swap chain image for directly
                                         // rendering, namely color attachment

    create_info.presentMode = mode;

    /*
        6.1 image sharing mode
        exclusive: ownership changed when an image is visited by multiple queue
        concurrent:
    */
    QueueFamilyIndices indices = findQueueFamilies(mPhysicalDevice, mSurface);
    uint32_t queueFamilyIndices[] = {indices.presentFamily.value(),
                                     indices.graphicsFamily.value()};
    if (indices.graphicsFamily == indices.presentFamily)
    {
        // the present, and graphics queue are the same, we can use the
        // exclusive mode because the ownership can never be changed
        create_info.imageSharingMode = VK_SHARING_MODE_EXCLUSIVE;
        create_info.queueFamilyIndexCount = 0;     // don't know
        create_info.pQueueFamilyIndices = nullptr; // don't know
    }
    else
    {
        // if they are different
        create_info.imageSharingMode = VK_SHARING_MODE_CONCURRENT;
        create_info.queueFamilyIndexCount = 2;                // don't know
        create_info.pQueueFamilyIndices = queueFamilyIndices; // don't know
    }

    // 6.2 set the pretrasform (90 degree clockwise rotation,etc)
    create_info.preTransform = details.capabilities.currentTransform;

    create_info.compositeAlpha =
        VK_COMPOSITE_ALPHA_OPAQUE_BIT_KHR; // blend the color with other window
    create_info.clipped = VK_TRUE;         // window, color obscured related...

    create_info.oldSwapchain =
        VK_NULL_HANDLE; // only one swapchain, so set it to false

    // 6.3 final create
    if (VK_SUCCESS !=
        vkCreateSwapchainKHR(mDevice, &create_info, nullptr, &mSwapChain))
    {
        SIM_ERROR("create swap chain failed");
    }
    // std::cout << "succ to create swap chain\n";
    uint32_t swapchain_image_count = 0;
    vkGetSwapchainImagesKHR(mDevice, mSwapChain, &swapchain_image_count,
                            nullptr);
    mSwapChainImages.resize(swapchain_image_count);
    vkGetSwapchainImagesKHR(mDevice, mSwapChain, &swapchain_image_count,
                            mSwapChainImages.data());
    mSwapChainImageFormat = format.format;
    mSwapChainExtent = extent;
    // SIM_INFO("swapchain is created successfully");
}

/**
 * \brief           Create image view
 */
VkImageView CreateImageView(VkDevice device, VkImage image, VkFormat format,
                            VkImageAspectFlags aspectFlags)
{
    VkImageViewCreateInfo viewInfo{};
    viewInfo.sType = VK_STRUCTURE_TYPE_IMAGE_VIEW_CREATE_INFO;
    viewInfo.image = image;
    viewInfo.viewType = VK_IMAGE_VIEW_TYPE_2D;
    viewInfo.format = format;
    viewInfo.subresourceRange.aspectMask = aspectFlags;
    viewInfo.subresourceRange.baseMipLevel = 0;
    viewInfo.subresourceRange.levelCount = 1;
    viewInfo.subresourceRange.baseArrayLayer = 0;
    viewInfo.subresourceRange.layerCount = 1;

    VkImageView imageView;
    if (vkCreateImageView(device, &viewInfo, nullptr, &imageView) != VK_SUCCESS)
    {
        throw std::runtime_error("failed to create texture image view!");
    }

    return imageView;
}

/**
 * \brief           create "view" for swap chain images
 *
 *          Frankly, we need to know how the SwapChainImages is accessed
 * */
void cDrawScene::CreateImageViews()
{
    mSwapChainImageViews.resize(mSwapChainImages.size());
    for (int i = 0; i < mSwapChainImageViews.size(); i++)
    {
        mSwapChainImageViews[i] =
            CreateImageView(mDevice, mSwapChainImages[i], mSwapChainImageFormat,
                            VK_IMAGE_ASPECT_COLOR_BIT);
    }
}

VkShaderModule cDrawScene::CreateShaderModule(const std::vector<char> &code)
{
    VkShaderModuleCreateInfo createInfo{};
    createInfo.sType = VK_STRUCTURE_TYPE_SHADER_MODULE_CREATE_INFO;
    createInfo.codeSize = code.size();
    createInfo.pCode = reinterpret_cast<const uint32_t *>(code.data());
    VkShaderModule module;
    // SIM_ASSERT(VK_SUCCESS ==);
    vkCreateShaderModule(mDevice, &createInfo, nullptr, &module);
    return module;
}

/**
 * \brief           Given physical device and tiling/feature requirements,
 * select the cnadidates VkFormat used for depth attachment.
 */
VkFormat findSupportedFormat(VkPhysicalDevice mPhyDevice,
                             const std::vector<VkFormat> &candidates,
                             VkImageTiling tiling,
                             VkFormatFeatureFlags features)
{
    for (const auto &format : candidates)
    {
        VkFormatProperties props;
        vkGetPhysicalDeviceFormatProperties(mPhyDevice, format, &props);

        // tiling feature
        if (tiling == VK_IMAGE_TILING_LINEAR &&
            (props.linearTilingFeatures & features) == features)
        {
            return format;
        }
        else if (tiling == VK_IMAGE_TILING_OPTIMAL &&
                 (props.optimalTilingFeatures & features) == features)
        {
            return format;
        }
    }
    exit(0);
}

VkFormat findDepthFormat(VkPhysicalDevice phy_device)
{
    /*
        For these candidates:
        VK_FORMAT_D32_SFLOAT: only depth buffer
        VK_FORMAT_D32_SFLOAT_S8_UINT: both depth and stencil buffer
        VK_FORMAT_D24_UNORM_S8_UINT: both depth and stencil buffer
    */
    // return findSupportedFormat(phy_device,
    //                            {VK_FORMAT_D32_SFLOAT,
    //                             VK_FORMAT_D32_SFLOAT_S8_UINT,
    //                             VK_FORMAT_D24_UNORM_S8_UINT},
    //                            VK_IMAGE_TILING_OPTIMAL,
    //                            VK_FORMAT_FEATURE_DEPTH_STENCIL_ATTACHMENT_BIT);
    return findSupportedFormat(phy_device, {VK_FORMAT_D32_SFLOAT},
                               VK_IMAGE_TILING_OPTIMAL,
                               VK_FORMAT_FEATURE_DEPTH_STENCIL_ATTACHMENT_BIT);
}

void cDrawScene::CreateRenderPass()
{
    // 1. describe the renderpass setting: an attachment is an image used in a
    // framebuffer.
    // the first color attachment (in MSAA)
    VkAttachmentDescription colorAttachment{};
    colorAttachment.format = mSwapChainImageFormat;
    colorAttachment.samples = mSampleCount; // no multisampling
    colorAttachment.loadOp =
        VK_ATTACHMENT_LOAD_OP_CLEAR; // clear the image in framebuffer, because
                                     // it is the first render pass
    colorAttachment.storeOp =
        VK_ATTACHMENT_STORE_OP_STORE; // store the result into the frame buffer
    colorAttachment.stencilLoadOp = VK_ATTACHMENT_LOAD_OP_DONT_CARE;   //
    colorAttachment.stencilStoreOp = VK_ATTACHMENT_STORE_OP_DONT_CARE; //
    colorAttachment.initialLayout =
        VK_IMAGE_LAYOUT_UNDEFINED; // which layout will the image have, before
                                   // render pass

    /*
        Note that for color image, instead of displaying directly, the
       multisampling image must be converted to the normal image. but for depth
       image, there is no such requirement.

       So here the finalLayout is changed
    */
    // colorAttachment.finalLayout =
    //     VK_IMAGE_LAYOUT_PRESENT_SRC_KHR; // which ... after render pass
    colorAttachment.finalLayout =
        VK_IMAGE_LAYOUT_COLOR_ATTACHMENT_OPTIMAL; // which ... after render pass

    // renference for attachment description, used for differnet subpasses
    VkAttachmentReference colorAttachmentRef{};
    colorAttachmentRef.attachment =
        0; // we have only one attachement, so it's 0
    colorAttachmentRef.layout = VK_IMAGE_LAYOUT_COLOR_ATTACHMENT_OPTIMAL;

    // 2. the second color resolved attachment (MSAA)
    VkAttachmentDescription colorAttachmentResolve{};
    colorAttachmentResolve.format = mSwapChainImageFormat;
    colorAttachmentResolve.samples = VK_SAMPLE_COUNT_1_BIT;
    colorAttachmentResolve.loadOp = VK_ATTACHMENT_LOAD_OP_DONT_CARE;
    colorAttachmentResolve.storeOp = VK_ATTACHMENT_STORE_OP_STORE;
    colorAttachmentResolve.stencilLoadOp = VK_ATTACHMENT_LOAD_OP_DONT_CARE;
    colorAttachmentResolve.stencilStoreOp = VK_ATTACHMENT_STORE_OP_DONT_CARE;
    colorAttachmentResolve.initialLayout = VK_IMAGE_LAYOUT_UNDEFINED;

    // if it is the last render pass, then the output layout must be presented

    if (mSceneIsTheLastRenderPass == true)
        colorAttachmentResolve.finalLayout = VK_IMAGE_LAYOUT_PRESENT_SRC_KHR;
    // if there is another render pass after it, we need to store the result
    // into the framebuffer
    else
        colorAttachmentResolve.finalLayout =
            VK_IMAGE_LAYOUT_COLOR_ATTACHMENT_OPTIMAL;

    VkAttachmentReference colorAttachmentResolveRef{};
    colorAttachmentResolveRef.attachment = 2;
    colorAttachmentResolveRef.layout = VK_IMAGE_LAYOUT_COLOR_ATTACHMENT_OPTIMAL;

    // add depth attachment
    VkAttachmentDescription depthAttachment{};
    depthAttachment.format = findDepthFormat(mPhysicalDevice);
    depthAttachment.samples = mSampleCount;
    depthAttachment.loadOp = VK_ATTACHMENT_LOAD_OP_CLEAR;
    depthAttachment.storeOp = VK_ATTACHMENT_STORE_OP_DONT_CARE;
    depthAttachment.stencilLoadOp = VK_ATTACHMENT_LOAD_OP_DONT_CARE;
    depthAttachment.stencilStoreOp = VK_ATTACHMENT_STORE_OP_DONT_CARE;
    depthAttachment.initialLayout = VK_IMAGE_LAYOUT_UNDEFINED;
    depthAttachment.finalLayout =
        VK_IMAGE_LAYOUT_DEPTH_STENCIL_ATTACHMENT_OPTIMAL;

    VkAttachmentReference depthAttachmentRef{};
    depthAttachmentRef.attachment = 1;
    depthAttachmentRef.layout =
        VK_IMAGE_LAYOUT_DEPTH_STENCIL_ATTACHMENT_OPTIMAL;

    // subpass
    VkSubpassDescription subpass{};
    subpass.pipelineBindPoint =
        VK_PIPELINE_BIND_POINT_GRAPHICS; // for display, not computing
    subpass.colorAttachmentCount = 1;
    subpass.pColorAttachments = &colorAttachmentRef;
    subpass.pDepthStencilAttachment = &depthAttachmentRef;
    subpass.pResolveAttachments = &colorAttachmentResolveRef;

    // where does this subpass come from? where will it go?
    VkSubpassDependency dependency{};
    dependency.srcSubpass = VK_SUBPASS_EXTERNAL;
    dependency.dstSubpass = 0;
    // support color and depth attachment
    dependency.srcStageMask = VK_PIPELINE_STAGE_COLOR_ATTACHMENT_OUTPUT_BIT |
                              VK_PIPELINE_STAGE_EARLY_FRAGMENT_TESTS_BIT;
    dependency.dstStageMask = VK_PIPELINE_STAGE_COLOR_ATTACHMENT_OUTPUT_BIT |
                              VK_PIPELINE_STAGE_EARLY_FRAGMENT_TESTS_BIT;
    dependency.dstAccessMask = VK_ACCESS_COLOR_ATTACHMENT_WRITE_BIT |
                               VK_ACCESS_DEPTH_STENCIL_ATTACHMENT_WRITE_BIT;

    // render pass
    std::array<VkAttachmentDescription, 3> attachments = {
        colorAttachment, depthAttachment, colorAttachmentResolve};

    VkRenderPassCreateInfo renderPassInfo{};
    renderPassInfo.sType = VK_STRUCTURE_TYPE_RENDER_PASS_CREATE_INFO;
    renderPassInfo.attachmentCount = attachments.size();
    renderPassInfo.pAttachments = attachments.data();
    renderPassInfo.subpassCount = 1;
    renderPassInfo.pSubpasses = &subpass;
    renderPassInfo.dependencyCount = 1;
    renderPassInfo.pDependencies = &dependency;
    SIM_ASSERT(vkCreateRenderPass(mDevice, &renderPassInfo, nullptr,
                                  &mRenderPass) == VK_SUCCESS);
}

void cDrawScene::CreateFrameBuffers()
{
    // framebuffers size = image views size
    mSwapChainFramebuffers.resize(mSwapChainImageViews.size());

    for (int i = 0; i < mSwapChainImageViews.size(); i++)
    {
        std::array<VkImageView, 3> attachments = {
            mColorImageViewMSAA, mDepthImageView, mSwapChainImageViews[i]};
        // VkImageView attachments[] = {mSwapChainImageViews[i]};
        VkFramebufferCreateInfo framebufferInfo{};
        framebufferInfo.sType = VK_STRUCTURE_TYPE_FRAMEBUFFER_CREATE_INFO;
        framebufferInfo.renderPass = mRenderPass;
        framebufferInfo.attachmentCount = attachments.size();
        framebufferInfo.pAttachments = attachments.data();
        framebufferInfo.width = mSwapChainExtent.width;
        framebufferInfo.height = mSwapChainExtent.height;
        framebufferInfo.layers = 1;
        SIM_ASSERT(VK_SUCCESS ==
                   vkCreateFramebuffer(mDevice, &framebufferInfo, nullptr,
                                       &(mSwapChainFramebuffers[i])));
        // SIM_ASSERT;
    }
}

void cDrawScene::CreateCommandPool()
{
    // VkCommandPool mCommandPool;
    auto queue_families = findQueueFamilies(mPhysicalDevice, mSurface);
    VkCommandPoolCreateInfo info{};
    info.sType = VK_STRUCTURE_TYPE_COMMAND_POOL_CREATE_INFO;
    info.queueFamilyIndex = queue_families.graphicsFamily.value();
    SIM_ASSERT(vkCreateCommandPool(mDevice, &info, nullptr, &mCommandPool) ==
               VK_SUCCESS);
    // SIM_INFO("Create Command Pool succ");
}

const int MAX_FRAMES_IN_FLIGHT = 2;
/**
 * \brief           Create (two) semaphores
 */
void cDrawScene::CreateSemaphores()
{
    mImageAvailableSemaphore.resize(MAX_FRAMES_IN_FLIGHT);
    mRenderFinishedSemaphore.resize(MAX_FRAMES_IN_FLIGHT);
    minFlightFences.resize(MAX_FRAMES_IN_FLIGHT);
    mImagesInFlight.resize(mSwapChainImages.size(), VK_NULL_HANDLE);

    VkSemaphoreCreateInfo semaphore_info{};
    semaphore_info.sType = VK_STRUCTURE_TYPE_SEMAPHORE_CREATE_INFO;

    VkFenceCreateInfo fence_info{};
    fence_info.sType = VK_STRUCTURE_TYPE_FENCE_CREATE_INFO;
    fence_info.flags = VK_FENCE_CREATE_SIGNALED_BIT;
    for (int i = 0; i < MAX_FRAMES_IN_FLIGHT; i++)
    {
        SIM_ASSERT(
            (VK_SUCCESS == vkCreateSemaphore(mDevice, &semaphore_info, nullptr,
                                             &mImageAvailableSemaphore[i])) &&
            (VK_SUCCESS == vkCreateSemaphore(mDevice, &semaphore_info, nullptr,
                                             &mRenderFinishedSemaphore[i])) &&
            (vkCreateFence(mDevice, &fence_info, nullptr,
                           &minFlightFences[i]) == VK_SUCCESS));
    }
    // SIM_INFO("Create semaphores succ");
}

/**
 * \brief           recreate the swap chain when the window is resized
 */
void cDrawScene::RecreateSwapChain()
{
    // wait
    vkDeviceWaitIdle(mDevice);

    // clean
    CleanSwapChain();

    // recreate
    CreateSwapChain();
    CreateImageViews();
    CreateRenderPass();
    CreateGraphicsPipeline("triangle", mTriangleGraphicsPipeline);
    CreateGraphicsPipeline("line", mLinesGraphicsPipeline);
    CreateCommandPool();
    CreateColorResources();
    CreateDepthResources();
    CreateFrameBuffers();
    CreateMVPUniformBuffer();
    CreateDescriptorPool();
    CreateDescriptorSets();
    CreateCommandBuffers();

    mImagesInFlight.resize(mSwapChainImages.size(), VK_NULL_HANDLE);
}

std::vector<tVkVertex> axes_vertices = {
    // X
    {{0.0f, 0.0f, 0.0f},
     {1.0f, 0.0f, 0.0f, 1.0f},
     {0.0f, 0.0f, 0.0f},
     {std::nanf(""), std::nanf("")}},
    {{100.0f, 0.0f, 0.0f},
     {1.0f, 0.0f, 0.0f, 1.0f},
     {0.0f, 0.0f, 0.0f},
     {std::nanf(""), std::nanf("")}},
    // Y
    {{0.0f, 0.0f, 0.0f},
     {0.0f, 1.0f, 0.0f, 1.0f},
     {0.0f, 0.0f, 0.0f},
     {std::nanf(""), std::nanf("")}},
    {{0.0f, 100.0f, 0.0f},
     {0.0f, 1.0f, 0.0f, 1.0f},
     {0.0f, 0.0f, 0.0f},
     {std::nanf(""), std::nanf("")}},
    // Z
    {{0.0f, 0.0f, 0.0f},
     {0.0f, 0.0f, 1.0f, 1.0f},
     {0.0f, 0.0f, 0.0f},
     {std::nanf(""), std::nanf("")}},
    {{0.0f, 0.0f, 100.0f},
     {0.0f, 0.0f, 1.0f, 1.0f},
     {0.0f, 0.0f, 0.0f},
     {std::nanf(""), std::nanf("")}},
};

void cDrawScene::CreateLineCommandBuffers(int i)
{
    vkCmdBindPipeline(mCommandBuffers[i], VK_PIPELINE_BIND_POINT_GRAPHICS,
                      mLinesGraphicsPipeline);

    // SIM_ERROR("VkBuffer mLineBuffer; vkDeviceMemory mLineBufferMemory; need
    // to be inited");
    VkBuffer lineBuffers[] = {mLineBuffer};
    VkDeviceSize offsets[] = {0};
    vkCmdBindVertexBuffers(mCommandBuffers[i], 0, 1, lineBuffers, offsets);

    vkCmdBindDescriptorSets(mCommandBuffers[i], VK_PIPELINE_BIND_POINT_GRAPHICS,
                            mPipelineLayout, 0, 1, &mDescriptorSets[i], 0,
                            nullptr);

    vkCmdDraw(mCommandBuffers[i], GetNumOfLineVertices(), 1, 0, 0);
}

void cDrawScene::CreateLineBuffer()
{
    VkDeviceSize buffer_size =
        sizeof(axes_vertices[0]) * GetNumOfLineVertices();
    if (buffer_size > SIM_VK_LINE_BUFFER_SIZE)
    {
        SIM_ERROR("current line buffer size {} > max size {}, please increase "
                  "SIM_VK_LINE_BUFFER_SIZE",
                  buffer_size, SIM_VK_LINE_BUFFER_SIZE);
        exit(1);
    }
    CreateBuffer(buffer_size, VK_BUFFER_USAGE_VERTEX_BUFFER_BIT,
                 VK_MEMORY_PROPERTY_HOST_VISIBLE_BIT |
                     VK_MEMORY_PROPERTY_HOST_COHERENT_BIT,
                 mLineBuffer, mLineBufferMemory);
}

void cDrawScene::UpdateLineBuffer(int idx)
{
    // update
    VkDeviceSize buffer_size =
        sizeof(axes_vertices[0]) * GetNumOfLineVertices();

    // 5. copy the vertex data to the buffer
    void *data = nullptr;
    // map the memory to "data" ptr;
    vkMapMemory(mDevice, mLineBufferMemory, 0, buffer_size, 0, &data);

    // write the data
    const tVectorXf &sim_edge_data = mSimScene->GetEdgesDrawBuffer();
    char *char_data = static_cast<char *>(data);
    memcpy(char_data, axes_vertices.data(),
           sizeof(axes_vertices[0]) * axes_vertices.size());
    memcpy(char_data + sizeof(axes_vertices[0]) * axes_vertices.size(),
           sim_edge_data.data(),
           sizeof(sim_edge_data[0]) * sim_edge_data.size());

    // unmap
    vkUnmapMemory(mDevice, mLineBufferMemory);
}

/**
 * \brief       judge whether a format contains the stencil buffer
 */
bool HasStencilComponent(VkFormat format)
{
    return format == VK_FORMAT_D32_SFLOAT_S8_UINT ||
           format == VK_FORMAT_D24_UNORM_S8_UINT;
}

uint32_t findMemoryType(VkPhysicalDevice physicalDevice, uint32_t typeFilter,
                        VkMemoryPropertyFlags properties)
{
    VkPhysicalDeviceMemoryProperties memProperties;
    vkGetPhysicalDeviceMemoryProperties(physicalDevice, &memProperties);

    for (uint32_t i = 0; i < memProperties.memoryTypeCount; i++)
    {
        if ((typeFilter & (1 << i)) &&
            (memProperties.memoryTypes[i].propertyFlags & properties) ==
                properties)
        {
            return i;
        }
    }

    throw std::runtime_error("failed to find suitable memory type!");
}

/**
 * \brief           create image
 */
void CreateImage(VkDevice device, VkPhysicalDevice phy_device, uint32_t width,
                 uint32_t height, VkSampleCountFlagBits num_samples,
                 VkFormat format, VkImageTiling tiling, VkImageUsageFlags usage,
                 VkMemoryPropertyFlags properties, VkImage &image,
                 VkDeviceMemory &imageMemory)
{
    VkImageCreateInfo imageInfo{};
    imageInfo.sType = VK_STRUCTURE_TYPE_IMAGE_CREATE_INFO;
    imageInfo.imageType = VK_IMAGE_TYPE_2D;
    imageInfo.extent.width = width;
    imageInfo.extent.height = height;
    imageInfo.extent.depth = 1;
    imageInfo.mipLevels = 1;
    imageInfo.arrayLayers = 1;
    imageInfo.format = format;
    imageInfo.tiling = tiling;
    imageInfo.initialLayout = VK_IMAGE_LAYOUT_UNDEFINED;
    imageInfo.usage = usage;
    imageInfo.samples = num_samples; // MSAA sampling
    imageInfo.sharingMode = VK_SHARING_MODE_EXCLUSIVE;

    if (vkCreateImage(device, &imageInfo, nullptr, &image) != VK_SUCCESS)
    {
        throw std::runtime_error("failed to create image!");
    }

    VkMemoryRequirements memRequirements;
    vkGetImageMemoryRequirements(device, image, &memRequirements);

    VkMemoryAllocateInfo allocInfo{};
    allocInfo.sType = VK_STRUCTURE_TYPE_MEMORY_ALLOCATE_INFO;
    allocInfo.allocationSize = memRequirements.size;
    allocInfo.memoryTypeIndex =
        findMemoryType(phy_device, memRequirements.memoryTypeBits, properties);

    if (vkAllocateMemory(device, &allocInfo, nullptr, &imageMemory) !=
        VK_SUCCESS)
    {
        throw std::runtime_error("failed to allocate image memory!");
    }

    vkBindImageMemory(device, image, imageMemory, 0);
}

extern VkCommandBuffer beginSingleTimeCommands(VkDevice device,
                                               VkCommandPool commandPool);
extern void endSingleTimeCommands(VkDevice device, VkCommandPool commandPool,
                                  VkQueue graphicsQueue,
                                  VkCommandBuffer commandBuffer);

void transitionImageLayout(VkDevice device, VkQueue graphicsQueue,
                           VkCommandPool commandPool, VkImage image,
                           VkFormat format, VkImageLayout oldLayout,
                           VkImageLayout newLayout)
{
    VkCommandBuffer commandBuffer =
        beginSingleTimeCommands(device, commandPool);

    VkImageMemoryBarrier barrier{};
    barrier.sType = VK_STRUCTURE_TYPE_IMAGE_MEMORY_BARRIER;
    barrier.oldLayout = oldLayout;
    barrier.newLayout = newLayout;
    barrier.srcQueueFamilyIndex = VK_QUEUE_FAMILY_IGNORED;
    barrier.dstQueueFamilyIndex = VK_QUEUE_FAMILY_IGNORED;
    barrier.image = image;
    barrier.subresourceRange.aspectMask = VK_IMAGE_ASPECT_COLOR_BIT;
    barrier.subresourceRange.baseMipLevel = 0;
    barrier.subresourceRange.levelCount = 1;
    barrier.subresourceRange.baseArrayLayer = 0;
    barrier.subresourceRange.layerCount = 1;

    VkPipelineStageFlags sourceStage;
    VkPipelineStageFlags destinationStage;

    if (newLayout == VK_IMAGE_LAYOUT_DEPTH_STENCIL_ATTACHMENT_OPTIMAL)
    {
        barrier.subresourceRange.aspectMask = VK_IMAGE_ASPECT_DEPTH_BIT;
    }
    else
    {
        barrier.subresourceRange.aspectMask = VK_IMAGE_ASPECT_COLOR_BIT;
    }

    // set up the source and destination mask
    if (oldLayout == VK_IMAGE_LAYOUT_UNDEFINED &&
        newLayout == VK_IMAGE_LAYOUT_TRANSFER_DST_OPTIMAL)
    {
        barrier.srcAccessMask = 0;
        barrier.dstAccessMask = VK_ACCESS_TRANSFER_WRITE_BIT;

        sourceStage = VK_PIPELINE_STAGE_TOP_OF_PIPE_BIT;
        destinationStage = VK_PIPELINE_STAGE_TRANSFER_BIT;
    }
    else if (oldLayout == VK_IMAGE_LAYOUT_TRANSFER_DST_OPTIMAL &&
             newLayout == VK_IMAGE_LAYOUT_SHADER_READ_ONLY_OPTIMAL)
    {
        barrier.srcAccessMask = VK_ACCESS_TRANSFER_WRITE_BIT;
        barrier.dstAccessMask = VK_ACCESS_SHADER_READ_BIT;

        sourceStage = VK_PIPELINE_STAGE_TRANSFER_BIT;
        destinationStage = VK_PIPELINE_STAGE_FRAGMENT_SHADER_BIT;
    }
    else if (oldLayout == VK_IMAGE_LAYOUT_UNDEFINED &&
             newLayout == VK_IMAGE_LAYOUT_DEPTH_STENCIL_ATTACHMENT_OPTIMAL)
    {
        barrier.srcAccessMask = 0;
        barrier.dstAccessMask = VK_ACCESS_DEPTH_STENCIL_ATTACHMENT_READ_BIT |
                                VK_ACCESS_DEPTH_STENCIL_ATTACHMENT_WRITE_BIT;

        sourceStage = VK_PIPELINE_STAGE_TOP_OF_PIPE_BIT;
        destinationStage = VK_PIPELINE_STAGE_EARLY_FRAGMENT_TESTS_BIT;
    }
    else
    {
        throw std::invalid_argument("unsupported layout transition!");
    }

    vkCmdPipelineBarrier(commandBuffer, sourceStage, destinationStage, 0, 0,
                         nullptr, 0, nullptr, 1, &barrier);

    endSingleTimeCommands(device, commandPool, graphicsQueue, commandBuffer);
}

/**
 * \brief           create depth resource
 */
void cDrawScene::CreateDepthResources()
{
    // 1. select a suitable depth format
    VkFormat depth_format = findDepthFormat(mPhysicalDevice);
    // std::cout << "[good] depth format = " << depth_format << std::endl;
    // create depth image
    CreateImage(
        mDevice, mPhysicalDevice, mSwapChainExtent.width,
        mSwapChainExtent.height, mSampleCount, depth_format,
        VK_IMAGE_TILING_OPTIMAL,
        // mSwapChainExtent.height, depth_format, VK_IMAGE_TILING_LINEAR,
        VK_IMAGE_USAGE_DEPTH_STENCIL_ATTACHMENT_BIT |
            VK_IMAGE_USAGE_TRANSFER_SRC_BIT,
        VK_MEMORY_PROPERTY_DEVICE_LOCAL_BIT, mDepthImage, mDepthImageMemory);

    // create image view
    mDepthImageView = CreateImageView(mDevice, mDepthImage, depth_format,
                                      VK_IMAGE_ASPECT_DEPTH_BIT);

    // cannot understand why we need to do transition
    transitionImageLayout(mDevice, mGraphicsQueue, mCommandPool, mDepthImage,
                          depth_format, VK_IMAGE_LAYOUT_UNDEFINED,
                          VK_IMAGE_LAYOUT_DEPTH_STENCIL_ATTACHMENT_OPTIMAL);
}

void cDrawScene::CreateColorResources()
{
    VkFormat color_format = mSwapChainImageFormat;
    CreateImage(mDevice, mPhysicalDevice, mSwapChainExtent.width,
                mSwapChainExtent.height, mSampleCount, color_format,
                VK_IMAGE_TILING_OPTIMAL,
                VK_IMAGE_USAGE_TRANSIENT_ATTACHMENT_BIT |
                    VK_IMAGE_USAGE_COLOR_ATTACHMENT_BIT,
                VK_MEMORY_PROPERTY_DEVICE_LOCAL_BIT, mColorImageMSAA,
                mColorImageMemoryMSAA);
    mColorImageViewMSAA = CreateImageView(
        mDevice, mColorImageMSAA, color_format, VK_IMAGE_ASPECT_COLOR_BIT);
}

/**
 * \brief               Create Image texture
 */
#define STB_IMAGE_IMPLEMENTATION
#include "utils/FileUtil.h"
#include "utils/stb_image.h"

void copyBufferToImage(VkDevice device, VkCommandPool commandPool,
                       VkQueue graphicsQueue, VkBuffer buffer, VkImage image,
                       uint32_t width, uint32_t height)
{
    VkCommandBuffer commandBuffer =
        beginSingleTimeCommands(device, commandPool);

    VkBufferImageCopy region{};
    region.bufferOffset = 0;
    region.bufferRowLength = 0;
    region.bufferImageHeight = 0;
    region.imageSubresource.aspectMask = VK_IMAGE_ASPECT_COLOR_BIT;
    region.imageSubresource.mipLevel = 0;
    region.imageSubresource.baseArrayLayer = 0;
    region.imageSubresource.layerCount = 1;
    region.imageOffset = {0, 0, 0};
    region.imageExtent = {width, height, 1};

    vkCmdCopyBufferToImage(commandBuffer, buffer, image,
                           VK_IMAGE_LAYOUT_TRANSFER_DST_OPTIMAL, 1, &region);

    endSingleTimeCommands(device, commandPool, graphicsQueue, commandBuffer);
}

void cDrawScene::CreateTextureImage()
{
    // 1. load the image from the file
    int tex_width, tex_height, tex_channels;
    // std::string mGroundPNGPath = "";
    SIM_ASSERT(cFileUtil::ExistsFile(mGroundPNGPath));
    stbi_uc *pixels = stbi_load(mGroundPNGPath.c_str(), &tex_width, &tex_height,
                                &tex_channels, STBI_rgb_alpha);
    SIM_ASSERT(tex_channels == 4);
    VkDeviceSize image_size = tex_width * tex_height * tex_channels;

    if (!pixels)
    {
        throw std::runtime_error("failed to load texture image!");
    }

    // 2. create the staging buffer, write the image to the staging buffer
    VkBuffer stagingBuffer;
    VkDeviceMemory stagingBufferMemory;
    CreateBuffer(image_size, VK_BUFFER_USAGE_TRANSFER_SRC_BIT,
                 VK_MEMORY_PROPERTY_HOST_VISIBLE_BIT |
                     VK_MEMORY_PROPERTY_HOST_COHERENT_BIT,
                 stagingBuffer, stagingBufferMemory);

    void *data = nullptr;
    vkMapMemory(mDevice, stagingBufferMemory, 0, image_size, 0, &data);
    memcpy(data, pixels, image_size);
    vkUnmapMemory(mDevice, stagingBufferMemory);

    stbi_image_free(pixels);

    // 3. create & allocate the texture image
    CreateImage(mDevice, mPhysicalDevice, tex_width, tex_height,
                VK_SAMPLE_COUNT_1_BIT, VK_FORMAT_R8G8B8A8_SRGB,
                VK_IMAGE_TILING_OPTIMAL,
                VK_IMAGE_USAGE_SAMPLED_BIT | VK_IMAGE_USAGE_TRANSFER_DST_BIT,
                VK_MEMORY_PROPERTY_DEVICE_LOCAL_BIT, mTextureImage,
                mTextureImageMemory);

    // 4. send the staging buffer to the image? how to do that?
    transitionImageLayout(mDevice, mGraphicsQueue, mCommandPool, mTextureImage,
                          VK_FORMAT_R8G8B8A8_SRGB, VK_IMAGE_LAYOUT_UNDEFINED,
                          VK_IMAGE_LAYOUT_TRANSFER_DST_OPTIMAL);
    copyBufferToImage(mDevice, mCommandPool, mGraphicsQueue, stagingBuffer,
                      mTextureImage, static_cast<uint32_t>(tex_width),
                      static_cast<uint32_t>(tex_height));
    transitionImageLayout(mDevice, mGraphicsQueue, mCommandPool, mTextureImage,
                          VK_FORMAT_R8G8B8A8_SRGB,
                          VK_IMAGE_LAYOUT_TRANSFER_DST_OPTIMAL,
                          VK_IMAGE_LAYOUT_SHADER_READ_ONLY_OPTIMAL);
    vkDestroyBuffer(mDevice, stagingBuffer, nullptr);
    vkFreeMemory(mDevice, stagingBufferMemory, nullptr);

    // SIM_INFO("Create texture image succ");
}

/**
 * \brief           create image view for texture
 */
void cDrawScene::CreateTextureImageView()
{
    mTextureImageView =
        CreateImageView(mDevice, mTextureImage, VK_FORMAT_R8G8B8A8_SRGB,
                        VK_IMAGE_ASPECT_COLOR_BIT);
    // SIM_INFO("create texture image view succ");
}

/**
 * \biref           create texture sampler and set this sampler
 */
void cDrawScene::CreateTextureSampler()
{
    VkSamplerCreateInfo samplerInfo{};
    samplerInfo.sType = VK_STRUCTURE_TYPE_SAMPLER_CREATE_INFO;
    samplerInfo.magFilter = VK_FILTER_LINEAR;
    samplerInfo.minFilter = VK_FILTER_LINEAR;
    samplerInfo.addressModeU = VK_SAMPLER_ADDRESS_MODE_REPEAT;
    samplerInfo.addressModeV = VK_SAMPLER_ADDRESS_MODE_REPEAT;
    samplerInfo.addressModeW = VK_SAMPLER_ADDRESS_MODE_REPEAT;
    samplerInfo.anisotropyEnable = VK_TRUE;
    VkPhysicalDeviceProperties properties{};
    vkGetPhysicalDeviceProperties(mPhysicalDevice, &properties);
    samplerInfo.maxAnisotropy = properties.limits.maxSamplerAnisotropy;
    samplerInfo.borderColor = VK_BORDER_COLOR_INT_OPAQUE_BLACK;
    samplerInfo.unnormalizedCoordinates = VK_FALSE;
    samplerInfo.compareEnable = VK_FALSE;
    samplerInfo.compareOp = VK_COMPARE_OP_ALWAYS;
    samplerInfo.mipmapMode = VK_SAMPLER_MIPMAP_MODE_LINEAR;
    samplerInfo.mipLodBias = 0.0f;
    samplerInfo.minLod = 0.0f;
    samplerInfo.maxLod = 0.0f;

    if (vkCreateSampler(mDevice, &samplerInfo, nullptr, &mTextureSampler) !=
        VK_SUCCESS)
    {
        throw std::runtime_error("failed to create texture sampler!");
    }
}
#include "utils/DefUtil.h"
int cDrawScene::GetNumOfLineVertices() const
{
    const tVectorXf &result = mSimScene->GetEdgesDrawBuffer();
    return axes_vertices.size() + result.size() / RENDERING_SIZE_PER_VERTICE;
}
