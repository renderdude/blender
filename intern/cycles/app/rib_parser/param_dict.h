#ifndef PARAM_DICT_H
#define PARAM_DICT_H

#include <map>
#include <unordered_map>

#include "util/types.h"

#include "error.h"
#include "parsed_parameter.h"

CCL_NAMESPACE_BEGIN

   enum class Parameter_Type
   {
      Boolean,
      Bxdf,
      Color,
      Integer,
      Normal,
      Parameter,
      Point2,
      Point3,
      Real,
      String,
      Texture,
      Vector2,
      Vector3
   };

   template < Parameter_Type PT >
   struct Parameter_Type_Traits
   {
   };

   class Parameter_Dictionary;
   using Mapped_Parameter_Dictionary =
       std::unordered_map< std::string, Parameter_Dictionary >;
   using Mapped_Vector_Dictionary =
       std::map< std::string, std::vector< Parameter_Dictionary > >;
   using Mapped_Attributes = std::unordered_map< std::string, Parsed_Parameter >;

   /** @brief Parameter_Dictionary.
    * @details
    */

   class Parameter_Dictionary
   {
    public:
      /// @name Initialization
      ///@{
      Parameter_Dictionary() = default;
      Parameter_Dictionary( Parsed_Parameter_Vector params );

      Parameter_Dictionary( Parsed_Parameter_Vector params0,
                            const Parsed_Parameter_Vector& params1 );
      ///@}

      /// @name Access
      ///@{
      std::string get_texture( const std::string& name ) const;
      const File_Loc* loc( const std::string& ) const;

      float get_one_float( const std::string& name, float def ) const;
      int get_one_int( const std::string& name, int def ) const;
      bool get_one_bool( const std::string& name, bool def ) const;

      float2 get_one_point2( const std::string& name, float2 def ) const;
      float2 get_one_vector2( const std::string& name, float2 def ) const;
      float3 get_one_color( const std::string& name, float3 def ) const;
      float3 get_one_point3( const std::string& name, float3 def ) const;
      float3 get_one_vector3( const std::string& name, float3 def ) const;
      float3 get_one_normal( const std::string& name, float3 def ) const;
      std::string get_one_string( const std::string& name,
                                  const std::string& def ) const;


      vector< float > get_float_array( const std::string& name ) const;
      vector< int > get_int_array( const std::string& name ) const;
      vector< uint8_t > get_bool_array( const std::string& name ) const;

      vector< float3 > get_color_array( const std::string& name ) const;
      vector< float2 > get_point2_array( const std::string& name ) const;
      vector< float2 > get_vector2_array( const std::string& name ) const;
      vector< float3 > get_point3_array( const std::string& name ) const;
      vector< float3 > get_vector3_array( const std::string& name ) const;
      vector< float3 > get_normal_array( const std::string& name ) const;
      vector< std::string > get_string_array( const std::string& name ) const;

      Parsed_Parameter const* get_parameter( const std::string& name ) const;
      Parsed_Parameter* get_parameter( const std::string& name );

      const Parsed_Parameter_Vector&
      get_parameter_vector() const
      {
         return params;
      }

      ///@}
      /// @name Status report
      ///@{
      void report_unused() const;
      ///@}
      /// @name Element change
      ///@{
      void rename_parameter( const std::string& before, const std::string& after );
      void rename_used_textures( const std::map< std::string, std::string >& m );
      ///@}
      /// @name Removal
      ///@{
      void free_parameters();
      void remove_float( const std::string& );
      void remove_int( const std::string& );
      void remove_bool( const std::string& );
      void remove_point2( const std::string& );
      void remove_vector2( const std::string& );
      void remove_point3( const std::string& );
      void remove_vector3( const std::string& );
      void remove_normal( const std::string& );
      void remove_string( const std::string& );
      void remove_texture( const std::string& );
      void remove_spectrum( const std::string& );
      void remove_bxdf( const std::string& );
      ///@}
      /// @name Conversion
      ///@{
      std::string to_string() const;
      ///@}
      /// @name Basic operations
      ///@{
      // RenderMan - Need to modify the parameters after they've been
      // created
      void
      push_back( Parsed_Parameter* param )
      {
         params.push_back( param );
      }
      ///@}
      /// @name Miscellaneous
      ///@{
      ///@}
      /// @name Obsolete
      ///@{
      ///@}
      /// @name Inapplicable
      ///@{
      ///@}

    protected:
      ///@{
      ///@}

    private:
      template < Parameter_Type PT >
      typename Parameter_Type_Traits< PT >::Return_Type lookup_single(
          const std::string& name,
          typename Parameter_Type_Traits< PT >::Return_Type default_value ) const;

      template < Parameter_Type PT >
      vector< typename Parameter_Type_Traits< PT >::Return_Type > lookup_array(
          const std::string& name ) const;

      template < typename Return_Type, typename G, typename C >
      vector< Return_Type > lookup_array( const std::string& name,
                                               Parameter_Type type,
                                               const char* typeName, int nPerItem,
                                               G getValues, C convert ) const;

      void remove( const std::string& name, const char* typeName );
      void check_parameter_types();
      static std::string to_parameter_definition( const Parsed_Parameter* p,
                                                  int indentCount );

      // Parameter_Dictionary Private Members
      Parsed_Parameter_Vector params;
      int _num_owned_params;
   };  // end of class Parameter_Dictionary

CCL_NAMESPACE_END

#endif  // PARAM_DICT_H
