#include "DrawSceneImGUI.h"
#include "SimScene.h"
#include "backends/imgui_impl_glfw.h"
#include "backends/imgui_impl_vulkan.h"
#include "imgui.h"
#include "utils/FileUtil.h"
#include "utils/LogUtil.h"
#include "utils/SysUtil.h"
#include "utils/TimeUtil.hpp"
#include <iostream>
#define GLFW_INCLUDE_VULKAN
#ifdef _WIN32
#define VK_USE_PLATFORM_WIN32_KHR
#include "GLFW/glfw3.h"
#include <format>
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
extern GLFWwindow *window;
extern QueueFamilyIndices findQueueFamilies(VkPhysicalDevice device,
                                            VkSurfaceKHR surface);
extern int GetMaxFramesInFlight();

// #ifdef IMGui_VULKAN_DEBUG_REPORT
static VKAPI_ATTR VkBool32 VKAPI_CALL
debug_report(VkDebugReportFlagsEXT flags, VkDebugReportObjectTypeEXT objectType,
             uint64_t object, size_t location, int32_t messageCode,
             const char *pLayerPrefix, const char *pMessage, void *pUserData)
{
    (void)flags;
    (void)object;
    (void)location;
    (void)messageCode;
    (void)pUserData;
    (void)pLayerPrefix; // Unused arguments
    fprintf(stderr,
            "[vulkan] Debug report from ObjectType: %i\nMessage: %s\n\n",
            objectType, pMessage);
    return VK_FALSE;
}
// #endif // IMGui_VULKAN_DEBUG_REPORT

cDrawSceneImGui::cDrawSceneImGui() { mSceneIsTheLastRenderPass = false; }

cDrawSceneImGui::~cDrawSceneImGui() {}

/**
 * \brief           add the debug report callback after the creation of an
 * instance
 */
void cDrawSceneImGui::CreateInstance()
{
    cDrawScene::CreateInstance();
    // auto vkCreateDebugReportCallbackEXT =
    //     (PFN_vkCreateDebugReportCallbackEXT)(vkGetInstanceProcAddr(
    //         mInstance, "vkCreateDebugReportCallbackEXT"));
    auto vkCreateDebugReportCallbackEXT = PFN_vkCreateDebugReportCallbackEXT(
        vkGetInstanceProcAddr(mInstance, "vkCreateDebugReportCallbackEXT"));
    SIM_ASSERT(vkCreateDebugReportCallbackEXT != nullptr);

    VkDebugReportCallbackCreateInfoEXT debug_report_ci = {};
    debug_report_ci.sType =
        VK_STRUCTURE_TYPE_DEBUG_REPORT_CALLBACK_CREATE_INFO_EXT;
    debug_report_ci.flags = VK_DEBUG_REPORT_ERROR_BIT_EXT |
                            VK_DEBUG_REPORT_WARNING_BIT_EXT |
                            VK_DEBUG_REPORT_PERFORMANCE_WARNING_BIT_EXT;
    debug_report_ci.pfnCallback = debug_report;
    debug_report_ci.pUserData = NULL;
    vkCreateDebugReportCallbackEXT(mInstance, &debug_report_ci, nullptr,
                                   &mDebugReport);
    // SIM_INFO("imgui create instance done");
    // exit(1);
}

static void check_vk_result(VkResult err)
{
    if (err == 0)
        return;
    fprintf(stderr, "[vulkan] Error: VkResult = %d\n", err);
    if (err < 0)
        abort();
}

void cDrawSceneImGui::CreateImGuiContext()
{
    IMGUI_CHECKVERSION();
    ImGui::CreateContext();
    ImGuiIO &io = ImGui::GetIO();
    (void)io;
    // io.ConfigFlags |= ImGuiConfigFlags_NavEnableKeyboard;     // Enable
    // Keyboard Controls io.ConfigFlags |= ImGuiConfigFlags_NavEnableGamepad; //
    // Enable Gamepad Controls

    // Setup Dear ImGui style
    ImGui::StyleColorsDark();
    // ImGui::StyleColorsClassic();

    const QueueFamilyIndices &queue_fami =
        findQueueFamilies(mPhysicalDevice, mSurface);

    // Setup Platform/Renderer backends
    ImGui_ImplGlfw_InitForVulkan(window, true);
    ImGui_ImplVulkan_InitInfo init_info = {};
    init_info.Instance = mInstance;
    init_info.PhysicalDevice = mPhysicalDevice;
    init_info.Device = mDevice;

    init_info.QueueFamily = queue_fami.graphicsFamily.value();
    init_info.Queue = mGraphicsQueue;
    init_info.PipelineCache = VK_NULL_HANDLE;
    init_info.DescriptorPool = mDescriptorPoolImGui;
    init_info.Allocator = nullptr;
    init_info.MinImageCount = 2;
    init_info.ImageCount = mSwapChainImages.size();
    // std::cout << "set imgui image count = " << init_info.ImageCount
    //           << std::endl;
    // exit(1);
    init_info.CheckVkResultFn = check_vk_result;
    ImGui_ImplVulkan_Init(&init_info, mRenderPassImGui);
    // std::cout << "init imgui done\n";
    // exit(1);
}
/**
 * \breif           Init vulkan
 */
void cDrawSceneImGui::InitVulkan()
{
    cDrawScene::InitVulkan();
    CreateRenderPassImGui();
    CreateCommandPoolImGui();
    CreateFrameBuffersImGui();
    CreateDescriptorPoolImGui();
    CreateCommandBuffersImGui();

    // configure for ImGui
    CreateImGuiContext();
    CreateFontImGui();
}

/**
 * \brief           ImGui needs it own complex descriptor
 */
void cDrawSceneImGui::CreateDescriptorPoolImGui()
{
    VkDescriptorPoolSize pool_sizes[] = {
        {VK_DESCRIPTOR_TYPE_SAMPLER, 1000},
        {VK_DESCRIPTOR_TYPE_COMBINED_IMAGE_SAMPLER, 1000},
        {VK_DESCRIPTOR_TYPE_SAMPLED_IMAGE, 1000},
        {VK_DESCRIPTOR_TYPE_STORAGE_IMAGE, 1000},
        {VK_DESCRIPTOR_TYPE_UNIFORM_TEXEL_BUFFER, 1000},
        {VK_DESCRIPTOR_TYPE_STORAGE_TEXEL_BUFFER, 1000},
        {VK_DESCRIPTOR_TYPE_UNIFORM_BUFFER, 1000},
        {VK_DESCRIPTOR_TYPE_STORAGE_BUFFER, 1000},
        {VK_DESCRIPTOR_TYPE_UNIFORM_BUFFER_DYNAMIC, 1000},
        {VK_DESCRIPTOR_TYPE_STORAGE_BUFFER_DYNAMIC, 1000},
        {VK_DESCRIPTOR_TYPE_INPUT_ATTACHMENT, 1000}};
    VkDescriptorPoolCreateInfo pool_info = {};
    pool_info.sType = VK_STRUCTURE_TYPE_DESCRIPTOR_POOL_CREATE_INFO;
    pool_info.flags = VK_DESCRIPTOR_POOL_CREATE_FREE_DESCRIPTOR_SET_BIT;
    pool_info.maxSets = 1000 * IM_ARRAYSIZE(pool_sizes);
    pool_info.poolSizeCount = (uint32_t)IM_ARRAYSIZE(pool_sizes);
    pool_info.pPoolSizes = pool_sizes;
    auto err = vkCreateDescriptorPool(mDevice, &pool_info, nullptr,
                                      &mDescriptorPoolImGui);
    check_vk_result(err);
    // std::cout << "create imgui descriptor pool done\n";
    // exit(1);
}

/**
 * \breif           The rendering process for ImGui must be seperated into
 * another render pass after render pass for the scene
 *
 *  We don't need MSAA in GUI rendering...
 */
void cDrawSceneImGui::CreateRenderPassImGui()
{
    VkAttachmentDescription colorAttachment{};
    colorAttachment.format = mSwapChainImageFormat;
    colorAttachment.samples = VK_SAMPLE_COUNT_1_BIT;
    colorAttachment.loadOp =
        VK_ATTACHMENT_LOAD_OP_LOAD; // load the image from buffer, do not clear
                                    // because now the scene has been rendered
                                    // into the framebuffer
    colorAttachment.storeOp =
        VK_ATTACHMENT_STORE_OP_STORE; // store the result into the
    colorAttachment.stencilLoadOp = VK_ATTACHMENT_LOAD_OP_DONT_CARE;
    colorAttachment.stencilStoreOp = VK_ATTACHMENT_STORE_OP_DONT_CARE;
    colorAttachment.initialLayout = VK_IMAGE_LAYOUT_COLOR_ATTACHMENT_OPTIMAL;

    // currently the ImGui render pass is the last one, so we directly present
    // it
    colorAttachment.finalLayout = VK_IMAGE_LAYOUT_PRESENT_SRC_KHR;

    // ---------------------- create attachment reference
    VkAttachmentReference color_attachment_ref{};
    color_attachment_ref.attachment = 0;
    color_attachment_ref.layout = VK_IMAGE_LAYOUT_COLOR_ATTACHMENT_OPTIMAL;

    // ----------------------- create subpass
    VkSubpassDescription subpass = {};
    subpass.pipelineBindPoint = VK_PIPELINE_BIND_POINT_GRAPHICS;
    subpass.colorAttachmentCount = 1;
    subpass.pColorAttachments = &color_attachment_ref;

    // ----------------------- connect the subpass with the previous one
    VkSubpassDependency dependency{};
    dependency.srcSubpass = VK_SUBPASS_EXTERNAL;
    dependency.dstSubpass = 0;
    dependency.srcStageMask = VK_PIPELINE_STAGE_COLOR_ATTACHMENT_OUTPUT_BIT;
    dependency.dstStageMask = VK_PIPELINE_STAGE_COLOR_ATTACHMENT_OUTPUT_BIT;
    dependency.srcAccessMask = 0;
    dependency.dstAccessMask = VK_ACCESS_COLOR_ATTACHMENT_WRITE_BIT;

    // ----------------------- create the final render pass
    VkRenderPassCreateInfo renderpassinfo{};
    renderpassinfo.sType = VK_STRUCTURE_TYPE_RENDER_PASS_CREATE_INFO;
    renderpassinfo.attachmentCount = 1;
    renderpassinfo.pAttachments = &colorAttachment;
    renderpassinfo.subpassCount = 1;
    renderpassinfo.pSubpasses = &subpass;
    renderpassinfo.dependencyCount = 1;
    renderpassinfo.pDependencies = &dependency;
    SIM_ASSERT(vkCreateRenderPass(mDevice, &renderpassinfo, nullptr,
                                  &mRenderPassImGui) == VK_SUCCESS);
}

VkCommandBuffer beginSingleTimeCommands(VkDevice device,
                                        VkCommandPool commandPool);
void endSingleTimeCommands(VkDevice device, VkCommandPool commandPool,
                           VkQueue graphicsQueue,
                           VkCommandBuffer commandBuffer);

/**
 * \brief        upload font to the GPU
 */
void cDrawSceneImGui::CreateFontImGui()
{
    VkCommandBuffer command_buffer =
        beginSingleTimeCommands(mDevice, mCommandPoolImGui);
    ImGui_ImplVulkan_CreateFontsTexture(command_buffer);
    endSingleTimeCommands(mDevice, mCommandPoolImGui, mGraphicsQueue,
                          command_buffer);
}

void cDrawSceneImGui::RecreateSwapChain()
{
    cDrawScene::RecreateSwapChain();

    // 1. all changed
    CreateRenderPassImGui();
    CreateCommandPoolImGui();

    CreateFrameBuffersImGui();
    CreateDescriptorPoolImGui();
    CreateCommandBuffersImGui();

    // configure for ImGui
    CreateImGuiContext();
    CreateFontImGui();
}
static void HelpMarker(const char *desc)
{
    ImGui::TextDisabled("(?)");
    if (ImGui::IsItemHovered())
    {
        ImGui::BeginTooltip();
        ImGui::PushTextWrapPos(ImGui::GetFontSize() * 35.0f);
        ImGui::TextUnformatted(desc);
        ImGui::PopTextWrapPos();
        ImGui::EndTooltip();
    }
}
extern int gWindowHeight, gWindowWidth;
void cDrawSceneImGui::Update(double dt)
{
    // update the GUI
    ImGui_ImplVulkan_NewFrame(); // empty
    ImGui_ImplGlfw_NewFrame();   // make sense
    ImGui::NewFrame();
    ImGui::GetIO().IniFilename = NULL; // disable ini saving

    const ImGuiViewport *main_viewport = ImGui::GetMainViewport();
    // std::cout << mSimScene->GetSceneType() << std::endl;
    if (mSimScene != nullptr)
    {
        switch (mSimScene->GetSceneType())
        {
        case SCENE_SIM:
        {
            UpdateSimImGui();
            break;
        }
        default:
        { // no gui
            std::cout << "Get scene type " << mSimScene->GetSceneType() << std::endl;
            break;
        }
        }
    }
    ImGui::Render();

    // update the simulation scene
    cDrawScene::Update(dt);
}

/**
 * \brief           imgui setting for acoustic simulation scene 
*/
void cDrawSceneImGui::UpdateSimImGui()
{
    ImVec2 init_window_size = ImVec2(400, 300);
    ImGui::SetNextWindowSize(init_window_size, ImGuiCond_FirstUseEver);

    ImGui::SetNextWindowPos(ImVec2(float(gWindowWidth) - init_window_size.x, 0),
                            ImGuiCond_FirstUseEver);

    ImGuiWindowFlags window_flags = 0;
    // window_flags |= ImGuiWindowFlags_NoMove;
    // window_flags |= ImGuiWindowFlags_NoResize;
    // bool open = false;
    // bool *p_open = &open;

    ImGui::Begin("config");
    // std::cout << "update acoustic imgui\n";
    // if (true == ImGui::Button("push_button"))
    // {
    //     std::cout << "[debug] button pushed\n";
    // }
    this->mSimScene->UpdateImGui();
    ImGui::End();
}
void cDrawSceneImGui::DrawFrame()
{
    uint32_t imageIndex;
    bool need_to_recreate_swapchain =
        FenceAndAcquireImageFromSwapchain(imageIndex);
    // updating the uniform buffer values
    // cTimeUtil::Begin("update_buffers");
    UpdateMVPUniformValue(imageIndex);
    UpdateVertexBufferSimObj(imageIndex);
    UpdateVertexBufferGround(imageIndex);
    UpdateLineBuffer(imageIndex);
    // cTimeUtil::End("update_buffers");

    // draw the GUI
    DrawImGui(imageIndex);

    // 2. submitting the command buffer
    std::vector<VkCommandBuffer> submitCommandBuffers = {
        mCommandBuffers[imageIndex], mCommandBufferImGui[imageIndex]};
    SubmitCmdsAndPresent(submitCommandBuffers, imageIndex);

    int MAX_FRAMES_IN_FLIGHT = GetMaxFramesInFlight();
    mCurFrame = (mCurFrame + 1) % MAX_FRAMES_IN_FLIGHT;

    if (need_to_recreate_swapchain == true)
    {
        mFrameBufferResized = false;
        RecreateSwapChain();
    }
}

/**
 * \brief           create a command pool for ImGui
 */
void cDrawSceneImGui::CreateCommandPoolImGui()
{
    VkCommandPoolCreateInfo info{};
    info.sType = VK_STRUCTURE_TYPE_COMMAND_POOL_CREATE_INFO;
    info.queueFamilyIndex =
        findQueueFamilies(mPhysicalDevice, mSurface).graphicsFamily.value();

    // this command pool has the ability to reset the command buffer
    info.flags = VK_COMMAND_POOL_CREATE_RESET_COMMAND_BUFFER_BIT;

    SIM_ASSERT(VK_SUCCESS == vkCreateCommandPool(mDevice, &info, nullptr,
                                                 &(mCommandPoolImGui)));
}

/**
 * \breif       create multiple command buffers for imgui
 */
void cDrawSceneImGui::CreateCommandBuffersImGui()
{
    uint32_t num_of_images = mSwapChainImages.size();
    mCommandBufferImGui.resize(num_of_images);
    VkCommandBufferAllocateInfo info{};
    info.sType = VK_STRUCTURE_TYPE_COMMAND_BUFFER_ALLOCATE_INFO;
    info.commandPool = mCommandPoolImGui;
    info.level = VK_COMMAND_BUFFER_LEVEL_PRIMARY;
    info.commandBufferCount = num_of_images;
    SIM_ASSERT(VK_SUCCESS == vkAllocateCommandBuffers(
                                 mDevice, &info, mCommandBufferImGui.data()));
}

/**
 * \brief      create the individual fraem buffers for imgui rendering
 */
void cDrawSceneImGui::CreateFrameBuffersImGui()
{
    // create image
    // create image view
    // create frame buffer

    uint32_t num_of_images = mSwapChainImages.size();
    mImGuiFrameBuffers.resize(num_of_images);
    VkImageView attachment[1];
    VkFramebufferCreateInfo info = {};
    info.sType = VK_STRUCTURE_TYPE_FRAMEBUFFER_CREATE_INFO;
    info.renderPass = this->mRenderPassImGui;
    info.attachmentCount = 1;
    info.pAttachments = attachment;
    info.width = mSwapChainExtent.width;
    info.height = mSwapChainExtent.height;
    info.layers = 1;

    for (int i = 0; i < mSwapChainImageViews.size(); i++)
    {
        attachment[0] = mSwapChainImageViews[i];
        SIM_ASSERT(vkCreateFramebuffer(mDevice, &info, nullptr,
                                       &(mImGuiFrameBuffers[i])) == VK_SUCCESS);
    }
}

void cDrawSceneImGui::DrawImGui(uint32_t imageIndex)
{

    auto err = vkResetCommandPool(mDevice, mCommandPoolImGui, 0);
    check_vk_result(err);
    VkCommandBufferBeginInfo begin_info = {};
    begin_info.sType = VK_STRUCTURE_TYPE_COMMAND_BUFFER_BEGIN_INFO;
    begin_info.flags |= VK_COMMAND_BUFFER_USAGE_ONE_TIME_SUBMIT_BIT;
    err = vkBeginCommandBuffer(mCommandBufferImGui[imageIndex], &begin_info);
    check_vk_result(err);

    VkClearValue clearValues;
    clearValues.color = {1.0f, 1.0f, 1.0f, 1.0f};

    VkRenderPassBeginInfo info = {};
    info.sType = VK_STRUCTURE_TYPE_RENDER_PASS_BEGIN_INFO;
    info.renderPass = mRenderPassImGui;
    info.framebuffer = mImGuiFrameBuffers[imageIndex];

    info.renderArea.extent.width = mSwapChainExtent.width;
    info.renderArea.extent.height = mSwapChainExtent.height;
    info.clearValueCount = 1;
    info.pClearValues = &clearValues;
    vkCmdBeginRenderPass(mCommandBufferImGui[imageIndex], &info,
                         VK_SUBPASS_CONTENTS_INLINE);

    ImGui_ImplVulkan_RenderDrawData(ImGui::GetDrawData(),
                                    mCommandBufferImGui[imageIndex]);
    // Submit command buffer
    vkCmdEndRenderPass(mCommandBufferImGui[imageIndex]);
    err = vkEndCommandBuffer(mCommandBufferImGui[imageIndex]);
    check_vk_result(err);
}

void cDrawSceneImGui::CleanSwapChain()
{
    cDrawScene::CleanSwapChain();

    // clean ImGui infos
    for (auto fb : mImGuiFrameBuffers)
    {
        vkDestroyFramebuffer(mDevice, fb, nullptr);
    }

    vkDestroyRenderPass(mDevice, mRenderPassImGui, nullptr);
    vkFreeCommandBuffers(mDevice, mCommandPoolImGui,
                         uint32_t(mCommandBufferImGui.size()),
                         mCommandBufferImGui.data());

    vkDestroyCommandPool(mDevice, mCommandPoolImGui, nullptr);

    // Resources to destroy when the program ends
    ImGui_ImplVulkan_Shutdown();
    ImGui_ImplGlfw_Shutdown();
    ImGui::DestroyContext();

    vkDestroyDescriptorPool(mDevice, mDescriptorPoolImGui, nullptr);
}

void cDrawSceneImGui::CursorMove(int xpos, int ypos)
{
    if (ImGui::GetIO().WantCaptureMouse == true)
        return;
    cDrawScene::CursorMove(xpos, ypos);
}
void cDrawSceneImGui::MouseButton(int button, int action, int mods)
{

    if (ImGui::GetIO().WantCaptureMouse == true)
        return;
    cDrawScene::MouseButton(button, action, mods);
}

/**
 * \brief       if Imgui want to control the scroll callback, we should stop to
 * callback it
 */
void cDrawSceneImGui::Scroll(double xoff, double yoff)
{
    if (ImGui::GetIO().WantCaptureMouse == true)
        return;
    else
    {
        cDrawScene::Scroll(xoff, yoff);
    }
}

/**
 * \brief           reset the scene imgui
 */
void cDrawSceneImGui::Reset()
{
    cDrawScene::Reset();
}

void cDrawSceneImGui::Init(const std::string &conf_path)
{
    cDrawScene::Init(conf_path);
}