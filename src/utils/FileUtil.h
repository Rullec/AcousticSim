#pragma once

#include "MathUtil.h"
#include "utils/LogUtil.h"
#include <fstream>
#include <string>

class cFileUtil
{
public:
    static FILE *OpenFile(const std::string &file_name, const char *mode);
    static FILE *OpenFile(const char *file_name, const char *mode);
    static void CloseFile(FILE *&f);
    static void ClearFile(const std::string &file_name);
    static void CreateFile(const std::string &file_name);
    static void DeleteFile(const char *file_name);
    static void ClearDir(const char *dir_name);
    static void DeleteDir(const char *dir_name);
    static bool CreateDir(const char *dir_name);
    static std::vector<std::string> ListDir(std::string);
    static std::string RemoveExtension(const std::string &filename);
    static void DeleteFile(const std::string &filename);
    static void RenameFile(const std::string &ori_name,
                           const std::string &des_name);
    static void CopyFile(const std::string &ori_name,
                         const std::string &des_name);
    static long int GetFileSize(const std::string &filename);
    static std::string GetExtension(const std::string &filename);
    static std::string GetFilename(const std::string &path);
    static std::string GetDir(const std::string &path);
    static std::string GetNoConflictDir(const std::string dir_base);
    static std::string GetNoConflictFile(const std::string file_base, const std::string &suffic);
    static void FilterFilesByExtension(std::vector<std::string> &files,
                                       const std::string &ext);
    static bool ExistsFile(const std::string &file_name);
    static bool ExistsDir(const std::string &dir_name);
    static bool ValidateFilePath(const std::string &file_name);
    static std::string ConcatFilename(const std::string &dir,
                                      const std::string &file);
    static void FindLine(std::ifstream &f_stream, int line);
    static std::string ReadTextFile(const std::string &path);

    // static bool ReadArray(FILE* f, const std::string& tag_beg, const
    // std::string& tag_end, std::vector<double>& out_buffer);
    static bool ReadTable(const std::string &filename,
                          std::vector<std::vector<double>> &out_buffer);
    static std::vector<std::string> ReadFileAllLines(const std::string &filename);
    static bool ReadMatrix(const std::string &filename,
                           Eigen::MatrixXd &out_mat);
    static bool WriteMatrix(const Eigen::MatrixXd &mat,
                            const std::string &out_filename);

    static bool AppendText(const std::string &str,
                           const std::string &out_filename);
    static std::string GenerateSerialFilename(const std::string &root, int id);
    static std::string GenerateRandomFilename(const std::string &root);

    // Multi Processes File Block
    static bool AddLock(const std::string &);
    static bool DeleteLock(const std::string &);

private:
    static std::string ReadTextFile(FILE *f);
    static std::string mFileLockDir;

    inline static const tLogger mLogger = cLogUtil::CreateLogger("FileUtil");
};
