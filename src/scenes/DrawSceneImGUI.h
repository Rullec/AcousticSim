#pragma once
#include "DrawScene.h"
#include "imgui.h"
class cDrawSceneImGui : public cDrawScene
{
public:
    explicit cDrawSceneImGui();
    virtual ~cDrawSceneImGui();
    virtual void Init(const std::string &conf_path) override;
    virtual void Update(double dt) override final;
    virtual void CursorMove(int xpos, int ypos) override;
    virtual void MouseButton(int button, int action, int mods) override;
    virtual void Scroll(double xoff, double yoff) override;
    virtual void Reset() override;

protected:
    virtual void InitVulkan() override;
    virtual void CreateInstance() override final;
    virtual void CreateImGuiContext();
    virtual void CreateDescriptorPoolImGui();
    virtual void CreateRenderPassImGui();
    virtual void CreateCommandPoolImGui();
    virtual void CreateCommandBuffersImGui();
    virtual void CreateFontImGui();
    virtual void RecreateSwapChain();
    virtual void CleanSwapChain() override;
    virtual void DrawFrame() override final;
    virtual void CreateFrameBuffersImGui();
    virtual void DrawImGui(uint32_t);
    virtual void UpdateSimImGui();
    VkDebugReportCallbackEXT mDebugReport;
    VkDescriptorPool mDescriptorPoolImGui;
    VkRenderPass mRenderPassImGui;
    VkCommandPool mCommandPoolImGui;
    

    std::vector<VkCommandBuffer> mCommandBufferImGui;
    std::vector<VkFramebuffer> mImGuiFrameBuffers;
    float mCurFPS;
};