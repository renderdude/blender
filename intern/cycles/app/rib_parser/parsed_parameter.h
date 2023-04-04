#ifndef PARSED_PARAMETER_H
#define PARSED_PARAMETER_H

#include <functional>
#include <memory>
#include <string>
#include <string_view>

#include "error.h"
#include "util/vector.h"

CCL_NAMESPACE_BEGIN
   /** @brief Parsed_Parameter.
    * @details
    */

   enum class Container_Type
   {
      Constant,
      FaceVarying,
      Reference,
      Uniform,
      Varying,
      Vertex
   };

   class Parsed_Parameter
   {
    public:
      /// @name Initialization
      ///@{
      Parsed_Parameter() = default;
      Parsed_Parameter( File_Loc loc ): loc( loc ) {}
      Parsed_Parameter( Parsed_Parameter const& pp );
      ///@}

      /// @name Access
      ///@{
      std::string type, name;
      int elem_per_item = 1;
      Container_Type storage = Container_Type::Constant;
      File_Loc loc;
      vector< float > floats;
      vector< int > ints;
      vector< std::string > strings;
      vector< uint8_t > bools;
      mutable bool looked_up                    = false;
      bool may_be_unused                        = false;
      ///@}
      /// @name Status Report
      ///@{
      bool
      floats_are( float value ) const
      {
         bool result = true;
         for ( auto it = floats.begin(); it != floats.end() && result; ++it )
            result &= ( *it == value );

         return result;
      }
      ///@}
      /// @name Conversion
      ///@{
      std::string to_string() const;
      ///@}
      /// @name Basic operations
      ///@{
      void add_float( float v );
      void add_int( int i );
      void add_string( std::string_view str );
      void add_bool( bool v );
      ///@}
   };  // end of class Parsed_Parameter

   using Parsed_Parameter_Vector = vector< Parsed_Parameter*>;

CCL_NAMESPACE_END

#endif  // PARSED_PARAMETER_H
