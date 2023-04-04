#ifndef PARSER_H
#define PARSER_H

#include <functional>
#include <map>
#include <memory>
#include <optional>
#include <set>
#include <string>
#include <string_view>
#include <unordered_map>

#include "error.h"
#include "param_dict.h"
#include "parsed_parameter.h"
#include "ri_api.h"
#include "scene_entities.h"

#include "util/span.h"
#include "util/transform.h"

CCL_NAMESPACE_BEGIN
using Point3f = float[3];
class Scene;

// Scene Parsing Declarations
void parse_files(Ri *target, std::vector<std::string> filenames);
void parse_string(Ri *target, std::string str);

// Token Definition
struct Token {
  Token() = default;
  Token(std::string_view token, File_Loc loc) : token(token), loc(loc)
  {
  }
  std::string to_string() const;
  std::string_view token;
  File_Loc loc;
};

// Tokenizer Definition
class Tokenizer {
 public:
  // Tokenizer Public Methods
  Tokenizer(std::string str,
            std::string filename,
            std::function<void(const char *, const File_Loc *)> error_callback);
#if defined(HAVE_MMAP) || defined(IS_WINDOWS)
  Tokenizer(void *ptr,
            size_t len,
            std::string filename,
            std::function<void(const char *, const File_Loc *)> error_callback);
#endif
  ~Tokenizer();

  static std::unique_ptr<Tokenizer> create_from_file(
      const std::string &filename,
      std::function<void(const char *, const File_Loc *)> error_callback);
  static std::unique_ptr<Tokenizer> create_from_string(
      std::string str, std::function<void(const char *, const File_Loc *)> error_callback);

  std::optional<Token> Next();

  // Just for parse().
  // TODO? Have a method to set this?
  File_Loc loc;

 private:
  void check_utf(const void *ptr, int len) const;
  // Tokenizer Private Methods
  int get_char()
  {
    if (pos == end)
      return EOF;
    int ch = *pos++;
    if (ch == '\n') {
      ++loc.line;
      loc.column = 0;
    }
    else
      ++loc.column;
    return ch;
  }
  void unget_char()
  {
    --pos;
    if (*pos == '\n')
      // Don't worry about the column; we'll be going to the start of
      // the next line again shortly...
      --loc.line;
  }

  // Tokenizer Private Members
  // This function is called if there is an error during lexing.
  std::function<void(const char *, const File_Loc *)> error_callback;

#if defined(HAVE_MMAP) || defined(IS_WINDOWS)
  // Scene files on disk are mapped into memory for lexing.  We need to
  // hold on to the starting pointer and total length so they can be
  // unmapped in the destructor.
  void *unmapPtr = nullptr;
  size_t unmapLength = 0;
#endif

  // If the input is stdin, then we copy everything until EOF into this
  // string and then start lexing.  This is a little wasteful (versus
  // tokenizing directly from stdin), but makes the implementation
  // simpler.
  std::string contents;

  // Pointers to the current position in the file and one past the end of
  // the file.
  const char *pos, *end;

  // If there are escaped characters in the string, we can't just return
  // a std::string_view into the mapped file. In that case, we handle the
  // escaped characters and return a std::string_view to sEscaped.  (And
  // thence, std::string_views from previous calls to Next() must be invalid
  // after a subsequent call, since we may reuse sEscaped.)
  std::string sEscaped;
};

CCL_NAMESPACE_END

#endif  // PARSER_H
