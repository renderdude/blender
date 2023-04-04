#ifndef ERROR_H
#define ERROR_H

#include <sstream>
#include <string>
#include <string_view>

CCL_NAMESPACE_BEGIN

struct File_Loc {
  File_Loc() = default;
  File_Loc(std::string_view filename) : filename(filename)
  {
  }

  std::string to_string() const
  {
    std::stringstream ss;
    ss << filename.data() << ":" << line << ":" << column;
    return ss.str();
  }

  std::string_view filename;
  int line = 1, column = 0;
};

inline void warning(const File_Loc *loc, std::string message)
{
  std::stringstream ss;
  ss << "Warning: " << loc->to_string() << ": " << message;
  fprintf(stderr, "%s\n", ss.str().c_str());
}

inline void error(const File_Loc *loc, std::string message)
{
  std::stringstream ss;
  ss << "Error: " << loc->to_string() << ": " << message;
  fprintf(stderr, "%s\n", ss.str().c_str());
}

inline void error_exit(const File_Loc *loc, std::string message)
{
  std::stringstream ss;
  ss << "Error: " << loc->to_string() << ": " << message;
  fprintf(stderr, "%s\n", ss.str().c_str());
  exit(EXIT_FAILURE);
}

CCL_NAMESPACE_END

#endif  // ERROR_H
