#include "app/cycles_xml.h"
#include "scene/camera.h"
#include "scene/scene.h"
#include "util/projection.h"
#include "util/transform.h"
#include <cstdlib>
#include <fcntl.h>
#include <sstream>
#include <sys/stat.h>

#include <fstream>
#include <functional>
#include <iostream>
#include <memory>
#include <string>
#include <utility>
#include <vector>

#include <boost/filesystem.hpp>
#include <double-conversion/double-conversion.h>
namespace bfs = boost::filesystem;

#include "util/log.h"

#include "error.h"
#include "parser.h"

CCL_NAMESPACE_BEGIN

static bfs::path search_directory;

void set_search_directory(std::string filename)
{
  bfs::path path(filename);
  if (!bfs::is_directory(path))
    path = path.parent_path();
  search_directory = path;
}

static std::string to_string_from_view(std::string_view s)
{
  return std::string(s.data(), s.size());
}

std::string Token::to_string() const
{
  return to_string_from_view(token);
}

// Tokenizer Implementation
static char decode_escaped(int ch, const File_Loc &loc)
{
  switch (ch) {
    case EOF:
      exit(-1);
    case 'b':
      return '\b';
    case 'f':
      return '\f';
    case 'n':
      return '\n';
    case 'r':
      return '\r';
    case 't':
      return '\t';
    case '\\':
      return '\\';
    case '\'':
      return '\'';
    case '\"':
      return '\"';
    default:
      exit(-1);
  }
  return 0;  // NOTREACHED
}

std::string read_file_contents(std::string filename)
{
#ifdef IS_WINDOWS
  std::ifstream ifs(WStringFromUTF8(filename).c_str(), std::ios::binary);
  if (!ifs)
    error_exit("%s: %s", filename, error_string());
  return std::string((std::istreambuf_iterator<char>(ifs)), (std::istreambuf_iterator<char>()));
#else
  int fd = open(filename.c_str(), O_RDONLY);
  if (fd == -1)
    exit(-1);

  struct stat stat;
  if (fstat(fd, &stat) != 0)
    exit(-1);

  std::string contents(stat.st_size, '\0');
  off_t chunk = 0;
  lseek(fd, 0, SEEK_SET);
  while (chunk < stat.st_size) {
    size_t readnow;
    readnow = read(fd, &contents[chunk], 1073741824);
    if (readnow < 0) 
      exit(-1);

    chunk = chunk + readnow;
  }
  close(fd);
  return contents;
#endif
}

static double_conversion::StringToDoubleConverter floatParser(
    double_conversion::StringToDoubleConverter::ALLOW_HEX,
    0. /* empty string value */,
    0. /* junk string value */,
    nullptr /* infinity symbol */,
    nullptr /* NaN symbol */);

std::unique_ptr<Tokenizer> Tokenizer::create_from_file(
    const std::string &filename,
    std::function<void(const char *, const File_Loc *)> error_callback)
{
  if (filename == "-") {
    // Handle stdin by slurping everything into a string.
    std::string str;
    int ch;
    while ((ch = getchar()) != EOF)
      str.push_back((char)ch);
    return std::make_unique<Tokenizer>(std::move(str), "<stdin>", std::move(error_callback));
  }

  // if ( has_extension( filename, ".gz" ) )
  //{
  //   std::string str = read_decompressed_file_contents( filename );
  //   return std::make_unique< Tokenizer >( std::move( str ), filename,
  //                                         std::move( error_callback ) );
  //}

#ifdef HAVE_MMAP
  int fd = open(filename.c_str(), O_RDONLY);
  if (fd == -1) {
    return nullptr;
  }

  struct stat stat;
  if (fstat(fd, &stat) != 0) {
    return nullptr;
  }

  size_t len = stat.st_size;
  if (len < 16 * 1024 * 1024) {
    close(fd);

    std::string str = read_file_contents(filename);
    return std::make_unique<Tokenizer>(std::move(str), filename, std::move(error_callback));
  }

  void *ptr = mmap(nullptr, len, PROT_READ, MAP_PRIVATE | MAP_NORESERVE, fd, 0);
  if (ptr == MAP_FAILED)

    if (close(fd) != 0) {
      return nullptr;
    }

  return std::make_unique<Tokenizer>(ptr, len, filename, std::move(error_callback));
#elif defined(IS_WINDOWS)
  auto errorReportLambda = [&error_callback, &filename]() -> std::unique_ptr<Tokenizer> {
    LPSTR messageBuffer = nullptr;
    FormatMessageA(FORMAT_MESSAGE_ALLOCATE_BUFFER | FORMAT_MESSAGE_FROM_SYSTEM |
                       FORMAT_MESSAGE_IGNORE_INSERTS,
                   NULL,
                   ::GetLastError(),
                   MAKELANGID(LANG_NEUTRAL, SUBLANG_DEFAULT),
                   (LPSTR)&messageBuffer,
                   0,
                   NULL);

    LocalFree(messageBuffer);
    return nullptr;
  };

  HANDLE fileHandle = CreateFileW(wstring_from_utf8(filename).c_str(),
                                  GENERIC_READ,
                                  FILE_SHARE_READ,
                                  0,
                                  OPEN_EXISTING,
                                  FILE_ATTRIBUTE_NORMAL,
                                  0);
  if (!fileHandle)
    return errorReportLambda();

  size_t len = GetFileSize(fileHandle, 0);

  HANDLE mapping = CreateFileMapping(fileHandle, 0, PAGE_READONLY, 0, 0, 0);
  CloseHandle(fileHandle);
  if (mapping == 0)
    return errorReportLambda();

  LPVOID ptr = MapViewOfFile(mapping, FILE_MAP_READ, 0, 0, 0);
  CloseHandle(mapping);
  if (!ptr) {
    return errorReportLambda();
  }

  return std::make_unique<Tokenizer>(ptr, len, filename, std::move(error_callback));
#else
  std::string str = read_file_contents(filename);
  return std::make_unique<Tokenizer>(std::move(str), filename, std::move(error_callback));
#endif
}

std::unique_ptr<Tokenizer> Tokenizer::create_from_string(
    std::string str, std::function<void(const char *, const File_Loc *)> error_callback)
{
  return std::make_unique<Tokenizer>(std::move(str), "<stdin>", std::move(error_callback));
}

Tokenizer::Tokenizer(std::string str,
                     std::string filename,
                     std::function<void(const char *, const File_Loc *)> error_callback)
    : error_callback(std::move(error_callback)), contents(std::move(str))
{
  loc = File_Loc(*new std::string(filename));
  pos = contents.data();
  end = pos + contents.size();
  check_utf(contents.data(), contents.size());
}

#if defined(HAVE_MMAP) || defined(IS_WINDOWS)
Tokenizer::Tokenizer(void *ptr,
                     size_t len,
                     std::string filename,
                     std::function<void(const char *, const File_Loc *)> error_callback)
    : error_callback(std::move(error_callback)), unmapPtr(ptr), unmapLength(len)
{
  // This is disgusting and leaks memory, but it ensures that the
  // filename in File_Locs returned by the Tokenizer remain valid even
  // after it has been destroyed.
  loc = File_Loc(*new std::string(filename));
  pos = (const char *)ptr;
  end = pos + len;
  check_utf(ptr, len);
}
#endif

Tokenizer::~Tokenizer()
{
#ifdef HAVE_MMAP
  if (unmapPtr && unmapLength > 0)
    if (munmap(unmapPtr, unmapLength) != 0)
      error_callback(string_printf("munmap: %s", error_string()).c_str(), nullptr);
#elif defined(IS_WINDOWS)
  if (unmapPtr && UnmapViewOfFile(unmapPtr) == 0)
    error_callback(string_printf("UnmapViewOfFile: %s", error_string()).c_str(), nullptr);
#endif
}

void Tokenizer::check_utf(const void *ptr, int len) const
{
  const unsigned char *c = (const unsigned char *)ptr;
  // https://en.wikipedia.org/wiki/Byte_order_mark
  if (len >= 2 && ((c[0] == 0xfe && c[1] == 0xff) || (c[0] == 0xff && c[1] == 0xfe)))
    error_callback(
        "File is encoded with UTF-16, which is not currently "
        "supported by pbrt (https://github.com/mmp/pbrt-v4/issues/136).",
        &loc);
}

std::optional<Token> Tokenizer::Next()
{
  while (true) {
    const char *tokenStart = pos;
    File_Loc startLoc = loc;

    int ch = get_char();
    if (ch == EOF)
      return {};
    else if (ch == ' ' || ch == '\n' || ch == '\t' || ch == '\r') {
      // nothing
    }
    else if (ch == '"') {
      // scan to closing quote
      bool haveEscaped = false;
      while ((ch = get_char()) != '"') {
        if (ch == EOF) {
          error_callback("premature EOF", &startLoc);
          return {};
        }
        else if (ch == '\n') {
          error_callback("unterminated string", &startLoc);
          return {};
        }
        else if (ch == '\\') {
          haveEscaped = true;
          // Grab the next character
          if ((ch = get_char()) == EOF) {
            error_callback("premature EOF", &startLoc);
            return {};
          }
        }
      }

      if (!haveEscaped)
        return Token({tokenStart, size_t(pos - tokenStart)}, startLoc);
      else {
        sEscaped.clear();
        for (const char *p = tokenStart; p < pos; ++p) {
          if (*p != '\\')
            sEscaped.push_back(*p);
          else {
            ++p;
            // CHECK_LT( p, pos );
            sEscaped.push_back(decode_escaped(*p, startLoc));
          }
        }
        return Token({sEscaped.data(), sEscaped.size()}, startLoc);
      }
    }
    else if (ch == '[' || ch == ']') {
      return Token({tokenStart, size_t(1)}, startLoc);
    }
    else if (ch == '#') {
      // comment: scan to EOL (or EOF)
      while ((ch = get_char()) != EOF) {
        if (ch == '\n' || ch == '\r') {
          unget_char();
          break;
        }
      }

      return Token({tokenStart, size_t(pos - tokenStart)}, startLoc);
    }
    else {
      // Regular statement or numeric token; scan until we hit a
      // space, opening quote, or bracket.
      while ((ch = get_char()) != EOF) {
        if (ch == ' ' || ch == '\n' || ch == '\t' || ch == '\r' || ch == '"' || ch == '[' ||
            ch == ']') {
          unget_char();
          break;
        }
      }
      return Token({tokenStart, size_t(pos - tokenStart)}, startLoc);
    }
  }
}

static int parse_int(const Token &t)
{
  bool negate = t.token[0] == '-';

  int index = 0;
  if (t.token[0] == '+' || t.token[0] == '-')
    ++index;

  int64_t value = 0;
  while (index < t.token.size()) {
    if (!(t.token[index] >= '0' && t.token[index] <= '9'))
      exit(-1);
    // error_exit( &t.loc, "\"%c\": expected a number", t.token[ index ] );
    value = 10 * value + (t.token[index] - '0');
    ++index;

    if (value > std::numeric_limits<int>::max())
      exit(-1);
    // error_exit(
    //    &t.loc,
    //    "Numeric value too large to represent as a 32-bit integer." );
    else if (value < std::numeric_limits<int>::lowest())
      std::cout << "Numeric value %d too low to represent as a 32-bit integer." << std::endl;
  }

  return negate ? -value : value;
}

static double parse_float(const Token &t)
{
  // Fast path for a single digit
  if (t.token.size() == 1) {
    if (!(t.token[0] >= '0' && t.token[0] <= '9'))
      exit(-1);
    // error_exit( &t.loc, "\"%c\": expected a number", t.token[ 0 ] );
    return t.token[0] - '0';
  }

  // Copy to a buffer so we can NUL-terminate it, as strto[idf]() expect.
  char buf[64];
  char *bufp = buf;
  std::unique_ptr<char[]> allocBuf;
  // CHECK_RARE( 1e-5, t.token.size() + 1 >= sizeof( buf ) );
  if (t.token.size() + 1 >= sizeof(buf)) {
    // This should be very unusual, but is necessary in case we get a
    // goofball number with lots of leading zeros, for example.
    allocBuf = std::make_unique<char[]>(t.token.size() + 1);
    bufp = allocBuf.get();
  }

  std::copy(t.token.begin(), t.token.end(), bufp);
  bufp[t.token.size()] = '\0';

  // Can we just use strtol?
  auto isInteger = [](std::string_view str) {
    for (char ch : str)
      if (!(ch >= '0' && ch <= '9'))
        return false;
    return true;
  };

  int length = 0;
  double val;
  if (isInteger(t.token)) {
    char *endptr;
    val = double(strtol(bufp, &endptr, 10));
    length = endptr - bufp;
  }
  else if (sizeof(float) == sizeof(float))
    val = (double)floatParser.StringToFloat(bufp, t.token.size(), &length);
  else
    val = floatParser.StringToDouble(bufp, t.token.size(), &length);

  if (length == 0)
    exit(-1);
  // error_exit( &t.loc, "%s: expected a number",
  //            to_string_from_view( t.token ) );

  return val;
}

inline bool is_quoted_string(std::string_view str)
{
  return str.size() >= 2 && str[0] == '"' && str.back() == '"';
}

static std::string_view dequote_string(const Token &t)
{
  if (!is_quoted_string(t.token))
    exit(-1);
  // error_exit( &t.loc, "\"%s\": expected quoted string",
  //            to_string_from_view( t.token ) );

  std::string_view str = t.token;
  str.remove_prefix(1);
  str.remove_suffix(1);
  return str;
}

constexpr int TokenOptional = 0;
constexpr int TokenRequired = 1;

template<typename Next, typename Unget>
static Parsed_Parameter_Vector parse_parameters(
    Next nextToken,
    Unget ungetToken,
    const std::function<void(const Token &token, const char *)> &error_callback)
{
  Parsed_Parameter_Vector parameterVector;

  while (true) {
    std::optional<Token> t = nextToken(TokenOptional);
    if (!t.has_value())
      return parameterVector;

    if (!is_quoted_string(t->token)) {
      ungetToken(*t);
      return parameterVector;
    }

    Parsed_Parameter *param = new Parsed_Parameter(t->loc);

    std::string_view decl = dequote_string(*t);

    auto skipSpace = [&decl](std::string_view::const_iterator iter) {
      while (iter != decl.end() && (*iter == ' ' || *iter == '\t'))
        ++iter;
      return iter;
    };
    // Skip to the next whitespace character (or the end of the string).
    auto skipToSpace = [&decl](std::string_view::const_iterator iter) {
      while (iter != decl.end() && *iter != ' ' && *iter != '\t')
        ++iter;
      return iter;
    };

    auto typeBegin = skipSpace(decl.begin());
    if (typeBegin == decl.end()) {
      std::cerr << &t->loc << "Parameter \"" << std::string(decl.begin(), decl.end())
                << "\" doesn't have a type declaration?!" << std::endl;
      exit(-1);
    }

    // Find end of type declaration
    auto typeEnd = skipToSpace(typeBegin);
    param->type.assign(typeBegin, typeEnd);

    // RenderMan
    bool was_storage = false;
    // Check to see if it's a container type
    if (param->type == "constant") {
      param->storage = Container_Type::Constant;
      was_storage = true;
    }
    else if (param->type == "facevarying") {
      param->storage = Container_Type::FaceVarying;
      was_storage = true;
    }
    else if (param->type == "reference") {
      param->storage = Container_Type::Reference;
      was_storage = true;
    }
    else if (param->type == "uniform") {
      param->storage = Container_Type::Uniform;
      was_storage = true;
    }
    else if (param->type == "varying") {
      param->storage = Container_Type::Varying;
      was_storage = true;
    }
    else if (param->type == "vertex") {
      param->storage = Container_Type::Vertex;
      was_storage = true;
    }

    if (was_storage) {
      typeBegin = skipSpace(typeEnd);
      if (typeBegin == decl.end()) {
        std::cerr << &t->loc << "Parameter \"" << std::string(decl.begin(), decl.end())
                  << "\" doesn't have a type declaration?!" << std::endl;
        exit(-1);
      }

      // Find end of type declaration
      typeEnd = skipToSpace(typeBegin);
      param->type.assign(typeBegin, typeEnd);
    }
    // See if it's an array
    std::vector<std::string> strings = split_string(param->type, '[');
    if (strings.size() > 1) {
      // Could extract the size and check we have that many values below
      param->type = strings[0];
      strings = split_string(strings[1], ']');
      param->elem_per_item = atoi(strings[0].c_str());
    }
    else {
      if (param->type == "point" || param->type == "normal")
        param->elem_per_item = 3;
      else if (param->type == "uv" || param->type == "st")
        param->elem_per_item = 2;
    }

    auto nameBegin = skipSpace(typeEnd);
    if (nameBegin == decl.end()) {
      std::cerr << &t->loc << " Unable to find parameter name from \""
                << std::string(decl.begin(), decl.end()) << "\"" << std::endl;
      exit(-1);
    }

    auto nameEnd = skipToSpace(nameBegin);
    param->name.assign(nameBegin, nameEnd);

    enum ValType { Unknown, String, Bool, Float, Int } valType = Unknown;

    if (param->type == "int")
      valType = Int;

    auto addVal = [&](const Token &t) {
      if (is_quoted_string(t.token)) {
        switch (valType) {
          case Unknown:
            valType = String;
            break;
          case String:
            break;
          case Float:
            error_callback(t, "expected floating-point value");
          case Int:
            error_callback(t, "expected integer value");
          case Bool:
            error_callback(t, "expected Boolean value");
        }

        param->add_string(dequote_string(t));
      }
      else if (t.token[0] == 't' && t.token == "true") {
        switch (valType) {
          case Unknown:
            valType = Bool;
            break;
          case String:
            error_callback(t, "expected string value");
          case Float:
            error_callback(t, "expected floating-point value");
          case Int:
            error_callback(t, "expected integer value");
          case Bool:
            break;
        }

        param->add_bool(true);
      }
      else if (t.token[0] == 'f' && t.token == "false") {
        switch (valType) {
          case Unknown:
            valType = Bool;
            break;
          case String:
            error_callback(t, "expected string value");
          case Float:
            error_callback(t, "expected floating-point value");
          case Int:
            error_callback(t, "expected integer value");
          case Bool:
            break;
        }

        param->add_bool(false);
      }
      else {
        switch (valType) {
          case Unknown:
            valType = Float;
            break;
          case String:
            error_callback(t, "expected string value");
          case Float:
            break;
          case Int:
            break;
          case Bool:
            error_callback(t, "expected Boolean value");
        }

        if (valType == Int)
          param->add_int(parse_int(t));
        else
          param->add_float(parse_float(t));
      }
    };

    Token val = *nextToken(TokenRequired);

    if (val.token == "[") {
      while (true) {
        val = *nextToken(TokenRequired);
        if (val.token == "]")
          break;
        addVal(val);
      }
    }
    else {
      addVal(val);
    }

    // RenderMan
    param->may_be_unused = true;

    parameterVector.push_back(param);
  }

  return parameterVector;
}

void parse(Ri *target, std::unique_ptr<Tokenizer> t)
{
  // LOG_VERBOSE( "Started parsing %s",
  //             std::string( t->loc.filename.begin(), t->loc.filename.end() ) );

  std::vector<std::unique_ptr<Tokenizer>> fileStack;
  fileStack.push_back(std::move(t));

  std::optional<Token> ungetToken;

  auto parse_error = [&](const char *msg, const File_Loc *loc) {
    std::cerr << loc << " " << msg << std::endl;
    exit(-2);
  };

  // nextToken is a little helper function that handles the file stack,
  // returning the next token from the current file until reaching EOF,
  // at which point it switches to the next file (if any).
  std::function<std::optional<Token>(int)> nextToken;
  nextToken = [&](int flags) -> std::optional<Token> {
    if (ungetToken.has_value())
      return std::exchange(ungetToken, {});

    if (fileStack.empty()) {
      if ((flags & TokenRequired) != 0) {
        std::cerr << "premature end of file" << std::endl;
        exit(-12);
      }
      return {};
    }

    std::optional<Token> tok = fileStack.back()->Next();

    if (!tok) {
      // We've reached EOF in the current file. Anything more to parse?
      // LOG_VERBOSE( "Finished parsing %s",
      //              std::string( fileStack.back()->loc.filename.begin(),
      //                           fileStack.back()->loc.filename.end() ) );
      fileStack.pop_back();
      return nextToken(flags);
    }
    else if (tok->token[0] == '#') {
      // Swallow comments, unless --format or --toply was given, in
      // which case they're printed to stdout.
      return nextToken(flags);
    }
    else
      // Regular token; success.
      return tok;
  };

  auto unget = [&](Token t) {
    // CHECK( !ungetToken.has_value() );
    ungetToken = t;
  };

  // Helper function for pbrt API entrypoints that take a single string
  // parameter and a ParameterVector (e.g. pbrtShape()).
  auto basic_param_list_entrypoint =
      [&](void (Ri::*apiFunc)(const std::string &, Parsed_Parameter_Vector, File_Loc),
          File_Loc loc) {
        Token t = *nextToken(TokenRequired);
        std::string_view dequoted = dequote_string(t);
        std::string n = to_string_from_view(dequoted);
        Parsed_Parameter_Vector parameterVector = parse_parameters(
            nextToken, unget, [&](const Token &t, const char *msg) {
              std::string token = to_string_from_view(t.token);
              std::string str = msg;
              parse_error(str.c_str(), &t.loc);
            });
        (target->*apiFunc)(n, std::move(parameterVector), loc);
      };

  // Helper function for RenderMan API entrypoints that take two string
  // parameters and a ParameterVector (e.g. Bxdf()).
  // using two_string_param_list_entrypoint = void (Basic_Scene_Builder::*)(const
  // std::string &, const std::string &, Parsed_Parameter_Vector, File_Loc);
  auto two_string_param_list_entrypoint =
      [&](void (Ri::*apiFunc)(
              const std::string &, const std::string &, Parsed_Parameter_Vector, File_Loc),
          File_Loc loc) {
        Token t = *nextToken(TokenRequired);
        std::string_view dequoted = dequote_string(t);
        std::string n = to_string_from_view(dequoted);
        t = *nextToken(TokenRequired);
        dequoted = dequote_string(t);
        std::string p = to_string_from_view(dequoted);
        Parsed_Parameter_Vector parameterVector = parse_parameters(
            nextToken, unget, [&](const Token &t, const char *msg) {
              std::string token = to_string_from_view(t.token);
              std::string str = msg;
              parse_error(str.c_str(), &t.loc);
            });
        (target->*apiFunc)(n, p, std::move(parameterVector), loc);
      };

  auto syntax_error = [&](const Token &t) {
    std::cerr << t.loc.to_string() << " Unknown directive: " << to_string_from_view(t.token)
              << std::endl;
    exit(-3);
    ;
  };

  auto parse_array_of_real = [&](std::vector<float> vals) {
    Token val = *nextToken(TokenRequired);

    if (val.token == "[") {
      while (true) {
        val = *nextToken(TokenRequired);
        if (val.token == "]")
          break;
        vals.push_back(parse_float(val));
      }
    }
    else
      syntax_error(val);
  };

  auto parse_array_of_int = [&](std::vector<int> &vals) {
    Token val = *nextToken(TokenRequired);

    if (val.token == "[") {
      while (true) {
        val = *nextToken(TokenRequired);
        if (val.token == "]")
          break;
        vals.push_back((int)(parse_int(val)));
      }
    }
    else
      syntax_error(val);
  };

  auto parse_array_of_string = [&](std::vector<std::string> &vals) {
    Token val = *nextToken(TokenRequired);

    if (val.token == "[") {
      while (true) {
        val = *nextToken(TokenRequired);
        if (val.token == "]")
          break;
        vals.push_back(to_string_from_view(dequote_string(val)));
      }
    }
    else
      syntax_error(val);
  };

  std::optional<Token> tok;

  while (true) {
    tok = nextToken(TokenOptional);
    if (!tok.has_value())
      break;

    switch (tok->token[0]) {
      case 'A':
        if (tok->token == "AttributeBegin")
          target->AttributeBegin(tok->loc);
        else if (tok->token == "AttributeEnd")
          target->AttributeEnd(tok->loc);
        else if (tok->token == "Attribute")
          basic_param_list_entrypoint(&Ri::Attribute, tok->loc);
        else if (tok->token == "AreaLightSource")
          basic_param_list_entrypoint(&Ri::AreaLightSource, tok->loc);
        // Ri API
        else if (tok->token == "Atmosphere")
          basic_param_list_entrypoint(&Ri::Atmosphere, tok->loc);
        else
          syntax_error(*tok);
        break;

      // Ri API
      case 'B':
        if (tok->token == "Basis") {
          Parsed_Parameter_Vector params = parse_parameters(
              nextToken, unget, [&](const Token &t, const char *msg) {
                std::string token = to_string_from_view(t.token);
                std::string str = msg;
                parse_error(str.c_str(), &t.loc);
              });
          target->Basis(std::move(params), tok->loc);
        }
        else if (tok->token == "Begin") {
          std::string_view n = dequote_string(*nextToken(TokenRequired));
          target->Begin(to_string_from_view(n), tok->loc);
        }
        else if (tok->token == "Bound") {
          bool bracketed = false;
          std::optional<Token> second = nextToken(TokenOptional);
          if (second.has_value()) {
            if (second->token == "[")
              bracketed = true;
            else
              syntax_error(*tok);
          }

          float m[6];
          for (int i = 0; i < 6; ++i)
            m[i] = parse_float(*nextToken(TokenRequired));
          if (bracketed)
            if (nextToken(TokenRequired)->token != "]")
              syntax_error(*tok);
          target->Bound(m, tok->loc);
        }
        else if (tok->token == "Bxdf")
          two_string_param_list_entrypoint(&Ri::Bxdf, tok->loc);
        else
          syntax_error(*tok);
        break;

      case 'C':
        if (tok->token == "ConcatTransform") {
          if (nextToken(TokenRequired)->token != "[")
            syntax_error(*tok);
          float m[16];
          for (int i = 0; i < 16; ++i)
            m[i] = parse_float(*nextToken(TokenRequired));
          if (nextToken(TokenRequired)->token != "]")
            syntax_error(*tok);
          target->ConcatTransform(m, tok->loc);
        }
        else if (tok->token == "CoordinateSystem") {
          std::string_view n = dequote_string(*nextToken(TokenRequired));
          target->CoordinateSystem(to_string_from_view(n), tok->loc);
        }
        else if (tok->token == "CoordSysTransform") {
          std::string_view n = dequote_string(*nextToken(TokenRequired));
          target->CoordSysTransform(to_string_from_view(n), tok->loc);
        }
        else if (tok->token == "Camera")
          basic_param_list_entrypoint(&Ri::camera, tok->loc);
        // Ri API
        else if (tok->token == "Clipping") {
          float m[2];
          for (int i = 0; i < 2; ++i)
            m[i] = parse_float(*nextToken(TokenRequired));
          target->Clipping(m[0], m[1], tok->loc);
        }
        else if (tok->token == "ClippingPlane") {
          float m[6];
          for (int i = 0; i < 6; ++i)
            m[i] = parse_float(*nextToken(TokenRequired));
          target->ClippingPlane(m[0], m[1], m[2], m[3], m[4], m[5], tok->loc);
        }
        else if (tok->token == "Color") {
          float m[3];
          for (int i = 0; i < 3; ++i)
            m[i] = parse_float(*nextToken(TokenRequired));
          target->Color(m[0], m[1], m[2], tok->loc);
        }
        else if (tok->token == "Cone") {
          bool bracketed = false;
          std::optional<Token> second = nextToken(TokenOptional);
          if (second.has_value()) {
            if (second->token == "[")
              bracketed = true;
            else
              unget(*second);
          }

          float m[3];
          for (int i = 0; i < 3; ++i)
            m[i] = parse_float(*nextToken(TokenRequired));

          if (bracketed)
            if (nextToken(TokenRequired)->token != "]")
              syntax_error(*tok);
          Parsed_Parameter_Vector params = parse_parameters(
              nextToken, unget, [&](const Token &t, const char *msg) {
                std::string token = to_string_from_view(t.token);
                std::string str = msg;
                parse_error(str.c_str(), &t.loc);
              });
          target->Cone(m[0], m[1], m[2], std::move(params), tok->loc);
        }
        else if (tok->token == "CropWindow") {
          float m[4];
          for (int i = 0; i < 4; ++i)
            m[i] = parse_float(*nextToken(TokenRequired));
          target->CropWindow(m[0], m[1], m[2], m[3], tok->loc);
        }
        else if (tok->token == "Curves") {
          std::string_view type = dequote_string(*nextToken(TokenRequired));

          std::vector<int> n_vertices;
          parse_array_of_int(n_vertices);

          std::string_view wrap = dequote_string(*nextToken(TokenRequired));
          Parsed_Parameter_Vector params = parse_parameters(
              nextToken, unget, [&](const Token &t, const char *msg) {
                std::string token = to_string_from_view(t.token);
                std::string str = msg;
                parse_error(str.c_str(), &t.loc);
              });
          target->Curves(to_string_from_view(type),
                         n_vertices,
                         to_string_from_view(wrap),
                         std::move(params),
                         tok->loc);
        }
        else if (tok->token == "Cylinder") {
          bool bracketed = false;
          std::optional<Token> second = nextToken(TokenOptional);

          if (second.has_value()) {
            if (second->token == "[")
              bracketed = true;
            else
              unget(*second);
          }

          float m[4];
          for (int i = 0; i < 4; ++i)
            m[i] = parse_float(*nextToken(TokenRequired));

          if (bracketed)
            if (nextToken(TokenRequired)->token != "]")
              syntax_error(*tok);
          Parsed_Parameter_Vector params = parse_parameters(
              nextToken, unget, [&](const Token &t, const char *msg) {
                std::string token = to_string_from_view(t.token);
                std::string str = msg;
                parse_error(str.c_str(), &t.loc);
              });
          target->Cylinder(m[0], m[1], m[2], m[3], std::move(params), tok->loc);
        }
        else
          syntax_error(*tok);
        break;

      // RI API
      case 'D':
        if (tok->token == "Declare") {
          std::string name = to_string_from_view(dequote_string(*nextToken(TokenRequired)));
          std::string decl = to_string_from_view(dequote_string(*nextToken(TokenRequired)));
          target->Declare(name, decl, tok->loc);
        }
        else if (tok->token == "DepthOfField") {
          float m[3];
          for (int i = 0; i < 3; ++i)
            m[i] = parse_float(*nextToken(TokenRequired));
          target->DepthOfField(m[0], m[1], m[2], tok->loc);
        }
        else if (tok->token == "Detail") {
          float m[6];
          for (int i = 0; i < 6; ++i)
            m[i] = parse_float(*nextToken(TokenRequired));
          target->Detail(m, tok->loc);
        }
        else if (tok->token == "DetailRange") {
          bool bracketed = false;
          std::optional<Token> second = nextToken(TokenOptional);
          if (second.has_value()) {
            if (second->token == "[")
              bracketed = true;
            else
              syntax_error(*tok);
          }

          float m[4];
          for (int i = 0; i < 4; ++i)
            m[i] = parse_float(*nextToken(TokenRequired));

          if (bracketed)
            if (nextToken(TokenRequired)->token != "]")
              syntax_error(*tok);
          target->DetailRange(m[0], m[1], m[2], m[3], tok->loc);
        }
        else if (tok->token == "Disk") {
          bool bracketed = false;
          std::optional<Token> second = nextToken(TokenOptional);
          if (second.has_value()) {
            if (second->token == "[")
              bracketed = true;
            else
              syntax_error(*tok);
          }

          float m[3];
          for (int i = 0; i < 3; ++i)
            m[i] = parse_float(*nextToken(TokenRequired));

          if (bracketed)
            if (nextToken(TokenRequired)->token != "]")
              syntax_error(*tok);
          Parsed_Parameter_Vector params = parse_parameters(
              nextToken, unget, [&](const Token &t, const char *msg) {
                std::string token = to_string_from_view(t.token);
                std::string str = msg;
                parse_error(str.c_str(), &t.loc);
              });
          target->Disk(m[0], m[1], m[2], std::move(params), tok->loc);
        }
        else if (tok->token == "Displace")
          two_string_param_list_entrypoint(&Ri::Displacement, tok->loc);
        else if (tok->token == "Display") {
          std::string_view name = dequote_string(*nextToken(TokenRequired));
          std::string_view type = dequote_string(*nextToken(TokenRequired));
          std::string_view mode = dequote_string(*nextToken(TokenRequired));
          Parsed_Parameter_Vector params = parse_parameters(
              nextToken, unget, [&](const Token &t, const char *msg) {
                std::string token = to_string_from_view(t.token);
                std::string str = msg;
                parse_error(str.c_str(), &t.loc);
              });
          target->Display(to_string_from_view(name),
                          to_string_from_view(type),
                          to_string_from_view(mode),
                          std::move(params),
                          tok->loc);
        }
        else if (tok->token == "DisplayChannel") {
          std::string_view channel = dequote_string(*nextToken(TokenRequired));
          Parsed_Parameter_Vector params = parse_parameters(
              nextToken, unget, [&](const Token &t, const char *msg) {
                std::string token = to_string_from_view(t.token);
                std::string str = msg;
                parse_error(str.c_str(), &t.loc);
              });
          target->DisplayChannel(to_string_from_view(channel), std::move(params), tok->loc);
        }
        else if (tok->token == "DisplayFilter") {
          std::string_view name = dequote_string(*nextToken(TokenRequired));
          std::string_view type = dequote_string(*nextToken(TokenRequired));
          Parsed_Parameter_Vector params = parse_parameters(
              nextToken, unget, [&](const Token &t, const char *msg) {
                std::string token = to_string_from_view(t.token);
                std::string str = msg;
                parse_error(str.c_str(), &t.loc);
              });
          target->DisplayFilter(
              to_string_from_view(name), to_string_from_view(type), std::move(params), tok->loc);
        }
        else
          syntax_error(*tok);
        break;

      // RI API
      case 'E':
        if (tok->token == "Else")
          target->Else(tok->loc);
        else if (tok->token == "ElseIf") {
          std::string_view expression = dequote_string(*nextToken(TokenRequired));
          target->ElseIf(to_string_from_view(expression), tok->loc);
        }
        else if (tok->token == "End")
          target->End(tok->loc);
        else if (tok->token == "Exposure") {
          float m[2];
          for (int i = 0; i < 2; ++i)
            m[i] = parse_float(*nextToken(TokenRequired));

          target->Exposure(m[0], m[1], tok->loc);
        }
        else if (tok->token == "Exterior") {
          basic_param_list_entrypoint(&Ri::Atmosphere, tok->loc);
        }
        else
          syntax_error(*tok);
        break;

      case 'F':
        if (tok->token == "Format") {
          int m[3];
          for (int i = 0; i < 3; ++i)
            m[i] = (int)parse_float(*nextToken(TokenRequired));

          target->Format(m[0], m[1], m[2], tok->loc);
        }
        else if (tok->token == "FrameAspectRatio")
          target->FrameAspectRatio(parse_float(*nextToken(TokenRequired)), tok->loc);
        else if (tok->token == "FrameBegin")
          target->FrameBegin((int)(parse_float(*nextToken(TokenRequired))), tok->loc);
        else if (tok->token == "FrameEnd")
          target->FrameEnd(tok->loc);
        else
          syntax_error(*tok);
        break;

      // RI API
      case 'G':
        if (tok->token == "GeneralPolygon") {
          std::vector<int> n_vertices;
          parse_array_of_int(n_vertices);

          Parsed_Parameter_Vector params = parse_parameters(
              nextToken, unget, [&](const Token &t, const char *msg) {
                std::string token = to_string_from_view(t.token);
                std::string str = msg;
                parse_error(str.c_str(), &t.loc);
              });
          target->GeneralPolygon(n_vertices, std::move(params), tok->loc);
        }
        else if (tok->token == "GeometricApproximation") {
          std::string_view type = dequote_string(*nextToken(TokenRequired));
          float m = parse_float(*nextToken(TokenRequired));

          target->GeometricApproximation(to_string_from_view(type), m, tok->loc);
        }
        else if (tok->token == "Geometry")
          basic_param_list_entrypoint(&Ri::Geometry, tok->loc);
        else
          syntax_error(*tok);
        break;

      // RI API
      case 'H':
        if (tok->token == "Hider") {
        }
        else if (tok->token == "Hyperboloid") {
        }
        else
          syntax_error(*tok);
        break;

      case 'I':
        if (tok->token == "Integrator")
          two_string_param_list_entrypoint(&Ri::Integrator, tok->loc);
        else if (tok->token == "Identity")
          target->Identity(tok->loc);
        // RI API
        else if (tok->token == "IfBegin") {
        }
        else if (tok->token == "IfEnd") {
        }
        else if (tok->token == "Illuminate") {
        }
        else if (tok->token == "Imager") {
        }
        else
          syntax_error(*tok);
        break;

      case 'L':
        if (tok->token == "LightSource")
          basic_param_list_entrypoint(&Ri::LightSource, tok->loc);
        // RI API
        else if (tok->token == "Light")
          two_string_param_list_entrypoint(&Ri::Light, tok->loc);
        else
          syntax_error(*tok);
        break;

      case 'M':
        if (tok->token == "MakeBrickMap") {
        }
        else if (tok->token == "MakeCubeFaceEnvironment") {
        }
        else if (tok->token == "MakeLatLongEnvironment") {
        }
        else if (tok->token == "MakeShadow") {
        }
        else if (tok->token == "MakeTexture") {
        }
        else if (tok->token == "Matte") {
        }
        else if (tok->token == "MotionBegin") {
        }
        else if (tok->token == "MotionEnd") {
        }
        else
          syntax_error(*tok);
        break;

      case 'N':
        if (tok->token == "NuPatch") {
        }
        else
          syntax_error(*tok);
        break;

      case 'O':
        if (tok->token == "ObjectBegin") {
          std::string_view n = dequote_string(*nextToken(TokenRequired));
          target->ObjectBegin(to_string_from_view(n), tok->loc);
        }
        else if (tok->token == "ObjectEnd")
          target->ObjectEnd(tok->loc);
        else if (tok->token == "ObjectInstance") {
          std::string_view n = dequote_string(*nextToken(TokenRequired));
          target->ObjectInstance(to_string_from_view(n), tok->loc);
        }
        else if (tok->token == "Option") {
          std::string_view n = dequote_string(*nextToken(TokenRequired));
          Parsed_Parameter_Vector params = parse_parameters(
              nextToken, unget, [&](const Token &t, const char *msg) {
                std::string token = to_string_from_view(t.token);
                std::string str = msg;
                parse_error(str.c_str(), &t.loc);
              });
          target->Option(to_string_from_view(n), std::move(params), tok->loc);
        }
        // RI API
        else if (tok->token == "Opacity") {
        }
        else if (tok->token == "Orientation") {
        }
        else
          syntax_error(*tok);
        break;

      case 'P':
        if (tok->token == "PixelFilter")
          basic_param_list_entrypoint(&Ri::PixelFilter, tok->loc);
        // RI API
        else if (tok->token == "Paraboloid") {
        }
        else if (tok->token == "Patch") {
        }
        else if (tok->token == "PatchMesh") {
        }
        else if (tok->token == "Pattern") {
          two_string_param_list_entrypoint(&Ri::Pattern, tok->loc);
        }
        else if (tok->token == "Perspective") {
        }
        else if (tok->token == "PixelSampleImager") {
        }
        else if (tok->token == "PixelVariance") {
        }
        else if (tok->token == "Points") {
        }
        else if (tok->token == "PointsGeneralPolygons") {
          std::vector<int> n_loops;
          parse_array_of_int(n_loops);

          std::vector<int> n_vertices;
          parse_array_of_int(n_vertices);

          std::vector<int> vertices;
          parse_array_of_int(vertices);

          Parsed_Parameter_Vector params = parse_parameters(
              nextToken, unget, [&](const Token &t, const char *msg) {
                std::string token = to_string_from_view(t.token);
                std::string str = msg;
                parse_error(str.c_str(), &t.loc);
              });
          target->PointsGeneralPolygons(
              n_loops, n_vertices, vertices, std::move(params), tok->loc);
        }
        else if (tok->token == "PointsPolygons") {
          std::vector<int> n_vertices;
          parse_array_of_int(n_vertices);

          std::vector<int> vertices;
          parse_array_of_int(vertices);

          Parsed_Parameter_Vector params = parse_parameters(
              nextToken, unget, [&](const Token &t, const char *msg) {
                std::string token = to_string_from_view(t.token);
                std::string str = msg;
                parse_error(str.c_str(), &t.loc);
              });
          target->PointsPolygons(n_vertices, vertices, std::move(params), tok->loc);
        }
        else if (tok->token == "Polygon") {
        }
        else if (tok->token == "ProcDelayedReadArchive") {
        }
        else if (tok->token == "Procedural") {
        }
        else if (tok->token == "Procedural2") {
        }
        else if (tok->token == "Projection")
          basic_param_list_entrypoint(&Ri::Projection, tok->loc);
        else
          syntax_error(*tok);
        break;

      // RI API
      case 'Q':
        if (tok->token == "Quantize") {
        }
        else
          syntax_error(*tok);
        break;

      case 'R':
        if (tok->token == "ReverseOrientation")
          target->ReverseOrientation(tok->loc);
        else if (tok->token == "Rotate") {
          float v[4];
          for (int i = 0; i < 4; ++i)
            v[i] = parse_float(*nextToken(TokenRequired));
          target->Rotate(v[0], v[1], v[2], v[3], tok->loc);
        }
        // RI API
        else if (tok->token == "ReadArchive") {
        }
        else if (tok->token == "RelativeDetail") {
        }
        else if (tok->token == "Resource") {
        }
        else if (tok->token == "ResourceBegin") {
        }
        else if (tok->token == "ResourceEnd") {
        }
        else
          syntax_error(*tok);
        break;

      case 'S':
        if (tok->token == "Scale") {
          float v[3];
          for (int i = 0; i < 3; ++i)
            v[i] = parse_float(*nextToken(TokenRequired));
          target->Scale(v[0], v[1], v[2], tok->loc);
        }
        // RI API
        else if (tok->token == "ScopedCoordinateSystem") {
          std::string_view n = dequote_string(*nextToken(TokenRequired));
          target->ScopedCoordinateSystem(to_string_from_view(n), tok->loc);
        }
        else if (tok->token == "ScreenWindow") {
        }
        else if (tok->token == "ShadingInterpolation") {
        }
        else if (tok->token == "ShadingRate") {
        }
        else if (tok->token == "Shutter") {
        }
        else if (tok->token == "Sides") {
        }
        else if (tok->token == "Shader") {
        }
        else if (tok->token == "Skew") {
        }
        else if (tok->token == "SolidBegin") {
        }
        else if (tok->token == "SolidEnd") {
        }
        else if (tok->token == "Sphere") {
          bool bracketed = false;
          std::optional<Token> second = nextToken(TokenOptional);

          if (second.has_value()) {
            if (second->token == "[")
              bracketed = true;
            else
              unget(*second);
          }

          float m[4];
          for (int i = 0; i < 4; ++i)
            m[i] = parse_float(*nextToken(TokenRequired));

          if (bracketed)
            if (nextToken(TokenRequired)->token != "]")
              syntax_error(*tok);
          Parsed_Parameter_Vector params = parse_parameters(
              nextToken, unget, [&](const Token &t, const char *msg) {
                std::string token = to_string_from_view(t.token);
                std::string str = msg;
                parse_error(str.c_str(), &t.loc);
              });
          target->Sphere(m[0], m[1], m[2], m[3], std::move(params), tok->loc);
        }
        else if (tok->token == "SubdivisionMesh") {
          std::string_view scheme = dequote_string(*nextToken(TokenRequired));
          int nfaces = parse_int(*nextToken(TokenRequired));

          std::vector<int> n_vertices;
          parse_array_of_int(n_vertices);

          std::vector<int> vertices;
          parse_array_of_int(vertices);

          std::vector<string> tags;
          parse_array_of_string(tags);

          std::vector<int> nargs;
          parse_array_of_int(nargs);

          std::vector<int> intargs;
          parse_array_of_int(intargs);

          std::vector<float> floatargs;
          parse_array_of_real(floatargs);

          Parsed_Parameter_Vector params = parse_parameters(
              nextToken, unget, [&](const Token &t, const char *msg) {
                std::string token = to_string_from_view(t.token);
                std::string str = msg;
                parse_error(str.c_str(), &t.loc);
              });
          target->SubdivisionMesh(to_string_from_view(scheme),
                                  nfaces,
                                  n_vertices,
                                  vertices,
                                  tags,
                                  nargs,
                                  intargs,
                                  floatargs,
                                  std::move(params),
                                  tok->loc);
        }
        else if (tok->token == "Surface") {
        }
        else if (tok->token == "System") {
        }
        else
          syntax_error(*tok);
        break;

      case 'T':
        if (tok->token == "TransformBegin") {
          {
            target->TransformBegin(tok->loc);
          }
        }
        else if (tok->token == "TransformEnd") {
          target->TransformEnd(tok->loc);
        }
        else if (tok->token == "Transform") {
          if (nextToken(TokenRequired)->token != "[")
            syntax_error(*tok);
          float m[16];
          for (int i = 0; i < 16; ++i)
            m[i] = parse_float(*nextToken(TokenRequired));
          if (nextToken(TokenRequired)->token != "]")
            syntax_error(*tok);
          target->transform(m, tok->loc);
        }
        else if (tok->token == "Translate") {
          float v[3];
          for (int i = 0; i < 3; ++i)
            v[i] = parse_float(*nextToken(TokenRequired));
          target->Translate(v[0], v[1], v[2], tok->loc);
        }
        else if (tok->token == "Texture") {
          std::string_view n = dequote_string(*nextToken(TokenRequired));
          std::string name = to_string_from_view(n);
          n = dequote_string(*nextToken(TokenRequired));
          std::string type = to_string_from_view(n);

          Token t = *nextToken(TokenRequired);
          std::string_view dequoted = dequote_string(t);
          std::string texName = to_string_from_view(dequoted);
          Parsed_Parameter_Vector params = parse_parameters(
              nextToken, unget, [&](const Token &t, const char *msg) {
                std::string token = to_string_from_view(t.token);
                std::string str = msg;
                parse_error(str.c_str(), &t.loc);
              });

          target->Texture(name, type, texName, std::move(params), tok->loc);
        }
        // RI API
        else if (tok->token == "TextureCoordinates") {
        }
        else if (tok->token == "Torus") {
        }
        else if (tok->token == "TransformPoints") {
        }
        else if (tok->token == "TrimCurve") {
        }
        else
          syntax_error(*tok);
        break;

      // RI API
      case 'v':
        if (tok->token == "version") {
          float version = parse_float(*nextToken(TokenRequired));
        }
        else
          syntax_error(*tok);
        break;

      case 'W':
        if (tok->token == "WorldBegin")
          target->WorldBegin(tok->loc);
        else if (tok->token == "WorldEnd")
          target->WorldEnd(tok->loc);
        else
          syntax_error(*tok);
        break;

      default:
        syntax_error(*tok);
    }
  }
}

void parse_files(Ri *target, std::vector<std::string> filenames)
{
  auto tok_error = [](const char *msg, const File_Loc *loc) {
    std::cerr << loc << " " << msg << std::endl;
    exit(-1);
  };

  // Process target description
  if (filenames.empty()) {
    // Parse target from standard input
    std::unique_ptr<Tokenizer> t = Tokenizer::create_from_file("-", tok_error);
    if (t)
      parse(target, std::move(t));
  }
  else {
    // Parse target from input files
    for (const std::string &fn : filenames) {
      if (fn != "-")
        set_search_directory(fn);

      std::unique_ptr<Tokenizer> t = Tokenizer::create_from_file(fn, tok_error);
      if (t)
        parse(target, std::move(t));
    }
  }

  target->end_of_files();
}

void parse_string(Ri *target, std::string str)
{
  auto tok_error = [](const char *msg, const File_Loc *loc) {
    std::cerr << loc << " " << msg << std::endl;
    exit(-1);
  };
  std::unique_ptr<Tokenizer> t = Tokenizer::create_from_string(std::move(str), tok_error);
  if (!t)
    return;
  parse(target, std::move(t));

  target->end_of_files();
}

CCL_NAMESPACE_END
