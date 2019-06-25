/**
 * @file
 * @brief Declares the VtkWriter which can write vtk files that can be read
 *        by ParaView
 * @author Raffael Casagrande
 * @date   2018-07-14 07:13:49
 * @copyright MIT License
 */

#ifndef __3e48c7b32a034cb3be3dbca884ff4f6c
#define __3e48c7b32a034cb3be3dbca884ff4f6c

#include <lf/mesh/mesh.h>
#include <lf/mesh/utils/utils.h>
#include <lf/uscalfe/uscalfe.h>
#include <Eigen/Eigen>
#include <boost/variant/variant.hpp>
#include <string>
#include <utility>
#include <vector>
#include "lf/mesh/utils/lambda_mesh_data_set.h"

namespace lf::io {

/**
 * @brief Representation of a VTK file (only relevant) features (advanced usage)
 *
 * This class can be written to a *.vtk file through the method WriteToFile().
 * The VtkWriter uses this class internally, it is not intended for direct
 * usage by the normal user but it can be useful if a special kind of
 * VtkFile should be written.
 */
class VtkFile {
 public:
  using size_type = unsigned int;

  /// Nested types
  //////////////////////////////////////////////////////////////////////////////
  enum class Format { ASCII, BINARY };

  enum class CellType {
    VTK_VERTEX = 1,
    VTK_POLY_VERTEX = 2,
    VTK_LINE = 3,
    VTK_POLY_LINE = 4,
    VTK_TRIANGLE = 5,
    VTK_TRIANGLE_STRIP = 6,
    VTK_POLYGON = 7,
    VTK_PIXEL = 8,
    VTK_QUAD = 9,
    VTK_TETRA = 10,
    VTK_VOXEL = 11,
    VTK_HEXAHEDRON = 12,
    VTK_WEDGE = 13,
    VTK_PYRAMID = 14,
    VTK_QUADRATIC_EDGE = 21,
    VTK_QUADRATIC_TRIANGLE = 22,
    VTK_QUADRATIC_QUAD = 23,
    VTK_QUADRATIC_TETRA = 24,
    VTK_QUADRATIC_HEXAHEDRON = 25,
    VTK_LAGRANGE_CURVE = 68,
    VTK_LAGRANGE_TRIANGLE = 69,
    VTK_LAGRANGE_QUADRILATERAL = 70,
    VTK_LAGRANGE_TETRAHEDRON = 71,
    VTK_LAGRANGE_HEXAHEDRON = 72,
    VTK_LAGRANGE_WEDGE = 73,
    VTK_LAGRANGE_PYRAMID = 74
  };

  class UnstructuredGrid {
   public:
    std::vector<Eigen::Vector3f> points;
    std::vector<std::vector<size_type>> cells;
    std::vector<CellType> cell_types;
  };

  template <class T>
  class FieldDataArray {
   public:
    std::string name;
    std::vector<T> data;

    FieldDataArray() = default;

    explicit FieldDataArray(std::string name, std::vector<T> data)
        : name(std::move(name)), data(std::move(data)) {}
  };

  /// Represents one set of attribute data (can be attached to points or cells)
  template <class T>
  class ScalarData {
   public:
    std::string name;
    std::vector<T> data;
    std::string lookup_table = "default";

    ScalarData() = default;

    explicit ScalarData(std::string data_name, std::vector<T> data,
                        std::string lookup_table = "default")
        : name(std::move(data_name)),
          data(std::move(data)),
          lookup_table(std::move(lookup_table)) {}
  };

  template <class T>
  class VectorData {
   public:
    std::string name;
    std::vector<Eigen::Matrix<T, 3, 1>> data;

    VectorData() = default;

    explicit VectorData(std::string name,
                        std::vector<Eigen::Matrix<T, 3, 1>> data)
        : name(std::move(name)), data(std::move(data)) {}
  };

  // WARNING: the type long is interepreted differently depending on the
  // platorm. I suggest to not use it at all.
  using Attributes = std::vector<boost::variant<
      ScalarData<char>, ScalarData<unsigned char>, ScalarData<short>,
      ScalarData<unsigned short>, ScalarData<int>, ScalarData<unsigned int>,
      ScalarData<long>, ScalarData<unsigned long>, ScalarData<float>,
      ScalarData<double>, VectorData<float>, VectorData<double>>>;

  using FieldData =
      std::vector<boost::variant<FieldDataArray<int>, FieldDataArray<float>,
                                 FieldDataArray<double>>>;

  // Actual members
  //////////////////////////////////////////////////////////////////////////////

  /// The Vtk Header, at most 256 characters, no new lines characters!
  std::string header;
  /// The format of the file.
  Format format = Format::ASCII;

  /// Describes the nodes + cells
  UnstructuredGrid unstructured_grid;

  FieldData field_data;

  // Data that is attached to points
  Attributes point_data;

  // Data that is attached to the cells of the mesh
  Attributes cell_data;
};

void WriteToFile(const VtkFile& vtk_file, const std::string& filename);

// clang-format off

/**
 * @brief Write a mesh along with mesh data into a vtk file.
 *
 *
 * ## Sample usage:
 * @snippet vtk_writer.cc usage
 *
 * ## Mesh
 * VtkWriter can work with any mesh implementing the `mesh::Mesh` interface.
 * By default, it will translate `codim=0` mesh entities into VTK cells.
 * If necessary, VtkWriter can also write out all `codim=1` mesh entities
 * (skeleton), see documentation of the Constructor.
 *
 *
 * ## Data
 * VtkWriter can attach different types of data to either cells (codim=0) or
 * points (codim=dimMesh) of the mesh (see the different overloads of
 * `WritePointData()` and `WriteCellData()`). Each dataset is translated into a
 * VTK AttributeSet which can then be visualized in ParaView.
 *
 * In addition to point- and cell-data VtkWriter can also write out global data
 * (in VTK terms this is called "field attribute sets"). This is useful e.g.
 * to write along the current time of the simulation and visualize it with a
 * label in ParaView.
 *
 *
 * @note The actual Vtk file is written when the VtkWriter is destructed!
 *       It is very important that this destructor is called.
 *
 * @note By default, VtkWriter outputs files in text-format. However, for
 * larger meshes you can switch to the more efficient binary mode (setBinary()).
 *
 * ## How to view the output in Paraview
 * Once a `*.vtk` file has been written with the VtkWriter, this file can be
 * opened in [Paraview](https://www.paraview.org/download/) as follows:
 *
 * 1. Open the file via `File->Open` (here we open the file `smiley.vtk` which
 * is produced by the example code in `examples/io/vtk_gmsh_demo`)
 * 2. Note that Paraview added the Vtk-File to the "Pipeline Browser". You can
 * activate the file and show it in the view by selecting the Vtk file in
 * the Pipeline Browser and clicking on "Apply" below:
 * ![](vtk_writer/paraview_apply.png)
 * 3. Paraview will now automatically select a set of data that was written with
 * `WritePointData()` or `WriteCellData()` and visualize it on the mesh. You
 * can select a different dataset or change the representation mode
 * (e.g. to show the underlying mesh):
 * ![](vtk_writer/paraview_data_selection.png)
 * 4. There are many ways to visualize vector datasets (i.e. if the MeshDataSet 
 * is `Eigen::Vector2d` or `Eigen::Vector3d` valued) in Paraview. If you
 * select a vector dataset, you can additionally select how the mesh should be
 * colored: Magnitude of the vector, or the x-, y- or z-component of the vector.
 * ![](vtk_writer/paraview_vector_dataset.png)
 * 5. If you want to visualize the vectors using arrows, you can use a
 * "Glyph" filter:
 *   1. Click on "Filters->Common->Glyph"
 *   2. Note that Paraview added the filter to the "Pipeline Browser"
 *   3. Select the new filter and adjust the properties in the "Properties Window"
 *   below:
 *   | Subsection   | Property Name | Description |
 *   | ------------ | ------------- | ----------- |
 *   | Glyph Source | Glyph Type    | `Arrow` is the most common choice, `2D Glyph` is also interesting
 *   | Active Attributes | Scalars  | The (scalar) dataset that is used to color the arrows. |
 *   | Active Attributes | Vector   | The vector dataset that defines the vector direction. If `Scalars` is a cell-based dataset, this must also be cell-based. If it `Scalars` is a point dataset, this must be a vector point dataset |
 *   | Scaling      | Scale Mode    | Set this to `vector` if the length of the arrows should be proportional to the magnitude of the vectors |
 *   | Scaling      | Scale Factor  | If vectors are too small/too big adjust this factor |
 *
 *   4. Click on "Apply"
 * 6. Finally, if you want to take a look at the raw data you can e.g. click on 
 * "Split Vertical" and then select "Spreadsheet View". You can then select 
 * rows and see to what cell/point they correspond in the mesh:
 * ![](vtk_writer/paraview_split_vertical.png)
 *
 * ## Higher order cells
 * Recent versions of Paraview [support higher order Lagrange cells](https://blog.kitware.com/modeling-arbitrary-order-lagrange-finite-elements-in-the-visualization-toolkit/)
 * of order up to 10. VtkWriter can write such higher order cells so that 
 * e.g. higher-order \ref Geometry "geometry mappings" or higher order BVP solutions can be visualized appropriately.
 * 
 * #### Example usage:
 * The following code will output a mesh of the unit circle (\f$ \{\vec{x} \in \mathbb{R}^2 | \norm{x}<1 \}\f$) consisting of second order quadrilaterals and the function \f$ \sin(\pi x)  \cos(\pi x)\f$ to vtk:
 * @snippet vtk_writer.cc highOrder
 * 
 * The two vtk files generated by this code can be visualized in Paraview:
 * ![](vtk_writer/2ndorder_comparison.png)
 * 
 * We can see that the first order mesh approximates the boundary of the circle with line segments and that the mesh function is approximated by bilinear basis functions.
 * In the present example this is not very satisfactory, but luckily second order cells remedy this situation.
 * 
 * Because it is rather hard to visualize second or higher order basis functions directly, Paraview interpolates such basis functions linearly. 
 * By default, Paraview approximates second order quadrilaterals by dividing them into four first order quadrilaterals. This is what is shown in the middle column above
 * where we can see that the boundary of the circle is still not very smooth. In order to get a better representation of the second order basis functions, we can instruct
 * Paraview to subdivide the cells even more with the tesselate filter (see below for instructions). This gives the desired result (on the right).
 * 
 * @note you need a recent version of Paraview to visualize higher-order cells. Version 5.6.0 is known to work.
 * 
 * @note By using higher order cells, not only the geometry is better approximated, the \emph Point datasets are also much better approximated. This is unfortunately not
 * true for cell datasets.
 * 
 * @note If you construct a VtkWriter with `order` greater than one, you cannot call WritePointData() methods anymore that accept a mesh::utils::MeshDataSet. 
 * If you want to write PointData, you need to call WritePointData(const std::string &name, const MESH_FUNCTION &mesh_function) which accepts a \ref mesh_function MeshFunction.
 * This is because a mesh::utils::MeshDataSet can only provide values at nodes of the mesh, but not at other (Lagrange) points that lie e.g. on the middle of an edge or 
 * inside of a cell. 
 * 
 * @note The current version of Paraview (5.6.0) has some issues with cells of [order 7 or higher and it just crashes](https://gitlab.kitware.com/paraview/paraview/issues/19022).
 *       Also there are [some problems with the visualization of higher order Lagrange curves](https://gitlab.kitware.com/paraview/paraview/issues/19023). This is only important
 *       when VtkWriter is constructed with `codim>0`.
 * 
 * #### How to apply the tesselate filter
 * As pointed out above, Paraview doesn't approximate higher-order Lagrange cells nicely by default. This can be resolved by applying the tesselate filter in Paraview:
 * 
 * 1) After you have loaded a VtkFile with `order>1` cells, click on "Filters->Search" in the menu bar. Then search for the "tesselate" filter.
 * 2) Make the following settings for the tesselate filter in the "properties" window:
 *   - "Chord error" determines the maximum geometric error. Set it e.g. to "0.001".
 *   - Just below the "Chord error" you can set the field error with a list of values. This determines the maximum visualization error of the datasets (e.g. the mesh function).
 *     It's not clear to the author (craffael) why we can specify multiple values here, but adding just one line with e.g. "0.0001" fits most purposes.   
  *    ![](vtk_writer/tesselate_settings.png)
 *  3) click on "Apply"
 */

// clang-format on
class VtkWriter {
 public:
  using dim_t = base::dim_t;

  VtkWriter(const VtkWriter&) = delete;
  VtkWriter(VtkWriter&&) = delete;
  VtkWriter& operator=(const VtkWriter&) = delete;
  VtkWriter& operator=(VtkWriter&&) = delete;

  /**
   * @brief Construct a new VtkWriter.
   * @param mesh The underlying mesh that should be written into the VtkFile.
   * @param filename The filename of the Vtk File
   * @param codim (Optional) the codimension of the cells, default is 0, i.e.
   *              the codim=0 entities of the mesh are written into the vtk
   *              files and data can be attached to these entities.
   *              If you set `codim=1`, the sekelton (only codim=1 entities)
   *              will be written into the vtk file and you can visualize data
   *              on the skeleton
   * @param order (Optional) the order of the cells written to Paraview. By
   *              default first order cells are used. Should be in [0,10], see
   *              also the documentation of the VtkWriter class for more info.
   */
  VtkWriter(std::shared_ptr<const mesh::Mesh> mesh, std::string filename,
            dim_t codim = 0, unsigned char order = 1);

  /**
   * @brief Determines whether the Vtk file is written in binary or ASCII mode
   *        (default).
   * @param binary  true if binary mode should be used, otherwise false
   */
  void setBinary(bool binary) {
    if (binary) {
      vtk_file_.format = VtkFile::Format::BINARY;
    } else {
      vtk_file_.format = VtkFile::Format::ASCII;
    }
  }

  /**
   * @brief Add a new `unsigned char` attribute dataset that attaches data to
   * the points/nodes of the mesh.
   * @param name The name of the attribute set.
   * @param mds The mesh dataset that that attaches the data to the points of
   *            the mesh.
   * @param undefined_value The value that should be written for a point to
   * which `mds` does not attach data (i.e. if `mds.DefinedOn() == false`)
   * @note This overload is only supported if VtkWriter has been instantiated
   * with `order=1`. For `order>1` use WritePointData(const std::string& <!--
   * --> name, const MESH_FUNCTION& mesh_function).
   */
  void WritePointData(const std::string& name,
                      const mesh::utils::MeshDataSet<unsigned char>& mds,
                      unsigned char undefined_value = 0);

  /**
   * @brief Add a new `char` attribute dataset that attaches data to
   * the points/nodes of the mesh.
   * @param name The name of the attribute set.
   * @param mds The mesh dataset that that attaches the data to the points of
   *            the mesh.
   * @param undefined_value The value that should be written for a point to
   * which `mds` does not attach data (i.e. if `mds.DefinedOn() == false`)
   * @note This overload is only supported if VtkWriter has been instantiated
   * with `order=1`. For `order>1` use WritePointData(const std::string& <!--
   * --> name, const MESH_FUNCTION& mesh_function).
   */
  void WritePointData(const std::string& name,
                      const mesh::utils::MeshDataSet<char>& mds,
                      char undefined_value = 0);

  /**
   * @brief Add a new `unsigned int` attribute dataset that attaches data to
   * the points/nodes of the mesh.
   * @param name The name of the attribute set.
   * @param mds The mesh dataset that that attaches the data to the points of
   *            the mesh.
   * @param undefined_value The value that should be written for a point to
   * which `mds` does not attach data (i.e. if `mds.DefinedOn() == false`)
   * @note This overload is only supported if VtkWriter has been instantiated
   * with `order=1`. For `order>1` use WritePointData(const std::string& <!--
   * --> name, const MESH_FUNCTION& mesh_function).
   */
  void WritePointData(const std::string& name,
                      const mesh::utils::MeshDataSet<unsigned int>& mds,
                      unsigned undefined_value = 0);

  /**
   * @brief Add a new `unsigned int` attribute dataset that attaches data to
   * the points/nodes of the mesh.
   * @param name The name of the attribute set.
   * @param mds The mesh dataset that that attaches the data to the points of
   *            the mesh.
   * @param undefined_value The value that should be written for a point to
   * which `mds` does not attach data (i.e. if `mds.DefinedOn() == false`)
   * @note This overload is only supported if VtkWriter has been instantiated
   * with `order=1`. For `order>1` use WritePointData(const std::string& <!--
   * --> name, const MESH_FUNCTION& mesh_function).
   */
  void WritePointData(const std::string& name,
                      const mesh::utils::MeshDataSet<int>& mds,
                      int undefined_value = 0);

  /**
   * @brief Add a new `float` attribute dataset that attaches data to
   * the points/nodes of the mesh.
   * @param name The name of the attribute set.
   * @param mds The mesh dataset that that attaches the data to the points of
   *            the mesh.
   * @param undefined_value The value that should be written for a point to
   * which `mds` does not attach data (i.e. if `mds.DefinedOn() == false`)
   * @note This overload is only supported if VtkWriter has been instantiated
   * with `order=1`. For `order>1` use WritePointData(const std::string& <!--
   * --> name, const MESH_FUNCTION& mesh_function).
   */
  void WritePointData(const std::string& name,
                      const mesh::utils::MeshDataSet<float>& mds,
                      float undefined_value = 0.f);

  /**
   * @brief Add a new `double` attribute dataset that attaches data to
   * the points/nodes of the mesh.
   * @param name The name of the attribute set.
   * @param mds The mesh dataset that that attaches the data to the points of
   *            the mesh.
   * @param undefined_value The value that should be written for a point to
   * which `mds` does not attach data (i.e. if `mds.DefinedOn() == false`)
   * @note This overload is only supported if VtkWriter has been instantiated
   * with `order=1`. For `order>1` use WritePointData(const std::string& <!--
   * --> name, const MESH_FUNCTION& mesh_function).
   */
  void WritePointData(const std::string& name,
                      const mesh::utils::MeshDataSet<double>& mds,
                      double undefined_value = 0.);

  /**
   * @brief Add a new vector attribute dataset that attaches vectors to
   * the points/nodes of the mesh.
   * @param name The name of the attribute set.
   * @param mds The mesh dataset that that attaches the data to the points of
   *            the mesh.
   * @param undefined_value The value that should be written for a point to
   * which `mds` does not attach data (i.e. if `mds.DefinedOn() == false`)
   * @note This overload is only supported if VtkWriter has been instantiated
   * with `order=1`. For `order>1` use WritePointData(const std::string& <!--
   * --> name, const MESH_FUNCTION& mesh_function).
   */
  void WritePointData(const std::string& name,
                      const mesh::utils::MeshDataSet<Eigen::Vector2d>& mds,
                      const Eigen::Vector2d& undefined_value = {0, 0});

  /**
   * @brief Add a new vector attribute dataset that attaches vectors to
   * the points/nodes of the mesh.
   * @param name The name of the attribute set.
   * @param mds The mesh dataset that that attaches the data to the points of
   *            the mesh.
   * @param undefined_value The value that should be written for a point to
   * which `mds` does not attach data (i.e. if `mds.DefinedOn() == false`)
   * @note This overload is only supported if VtkWriter has been instantiated
   * with `order=1`. For `order>1` use WritePointData(const std::string& <!--
   * --> name, const MESH_FUNCTION& mesh_function).
   */
  void WritePointData(const std::string& name,
                      const mesh::utils::MeshDataSet<Eigen::Vector2f>& mds,
                      const Eigen::Vector2f& undefined_value = {0, 0});

  /**
   * @brief Add a new vector attribute dataset that attaches vectors to
   * the points/nodes of the mesh.
   * @param name The name of the attribute set.
   * @param mds The mesh dataset that that attaches the data to the points of
   *            the mesh.
   * @param undefined_value The value that should be written for a point to
   * which `mds` does not attach data (i.e. if `mds.DefinedOn() == false`)
   * @note This overload is only supported if VtkWriter has been instantiated
   * with `order=1`. For `order>1` use WritePointData(const std::string& <!--
   * --> name, const MESH_FUNCTION& mesh_function).
   */
  void WritePointData(const std::string& name,
                      const mesh::utils::MeshDataSet<Eigen::Vector3d>& mds,
                      const Eigen::Vector3d& undefined_value = {0, 0, 0});

  /**
   * @brief Add a new vector attribute dataset that attaches vectors to
   * the points/nodes of the mesh.
   * @param name The name of the attribute set.
   * @param mds The mesh dataset that that attaches the data to the points
   * of the mesh.
   * @param undefined_value The value that should be written for a point to
   * which `mds` does not attach data (i.e. if `mds.DefinedOn() == false`)
   * @note This overload is only supported if VtkWriter has been instantiated
   * with `order=1`. For `order>1` use WritePointData(const std::string& <!--
   * --> name, const MESH_FUNCTION& mesh_function).
   */
  void WritePointData(const std::string& name,
                      const mesh::utils::MeshDataSet<Eigen::Vector3f>& mds,
                      const Eigen::Vector3f& undefined_value = {0, 0, 0});

  /**
   * @brief Add a new vector attribute dataset that attaches vectors to
   * the points/nodes of the mesh.
   * @param name The name of the attribute set.
   * @param mds The mesh dataset that that attaches the data to the points of
   *            the mesh.
   * @param undefined_value The value that should be written for a point to
   * which `mds` does not attach data (i.e. if `mds.DefinedOn() == false`)
   * @note This version accepts in principle arbitrary size double Vectors,
   * however only the first three components are visualized!
   * @note This overload is only supported if VtkWriter has been instantiated
   * with `order=1`. For `order>1` use WritePointData(const std::string& <!--
   * --> name, const MESH_FUNCTION& mesh_function).
   */
  void WritePointData(
      const std::string& name,
      const mesh::utils::MeshDataSet<Eigen::VectorXd>& mds,
      const Eigen::VectorXd& undefined_value = Eigen::Vector3d(0, 0, 0));

  /**
   * @brief Add a new vector attribute dataset that attaches vectors to
   * the points/nodes of the mesh.
   * @param name The name of the attribute set.
   * @param mds The mesh dataset that that attaches the data to the points of
   *            the mesh.
   * @param undefined_value The value that should be written for a point to
   * which `mds` does not attach data (i.e. if `mds.DefinedOn() == false`)
   * @note This version accepts in principle arbitrary size float Vectors,
   * however only the first three components are visualized!
   * @note This overload is only supported if VtkWriter has been instantiated
   * with `order=1`. For `order>1` use WritePointData(const std::string& <!--
   * --> name, const MESH_FUNCTION& mesh_function).
   */
  void WritePointData(
      const std::string& name,
      const mesh::utils::MeshDataSet<Eigen::VectorXf>& mds,
      const Eigen::VectorXf& undefined_value = Eigen::Vector3f(0, 0, 0));

  /**
   * @brief Sample a \ref mesh_function "MeshFunction" at points of the mesh and
   * possible edges/interior (for `order>1`).
   *
   * @tparam MESH_FUNCTION An object fulfilling the \ref mesh_function
   * concept. The \ref uscalfe::MeshFunctionReturnType should be one of
   * - unsigned char
   * - char
   * - unsigned int
   * - int
   * - float
   * - double
   * - Eigen::Vector2d
   * - Eigen::Vector2f
   * - Eigen::Vector3d
   * - Eigen::Vector3f
   * @param name The name of the dataset, cannot contain spaces!
   * @param mesh_function The \ref mesh_function "Mesh Function" to be sampled.
   *
   * @note this function will evaluate the MeshFunction on the points
   * of the mesh, i.e. at entities with codim=dimMesh. If the VtkWriter has been
   * instantiated with `order>1`, the mesh function will also be evaluated on
   * edges and possibly in the interior of the elements. Some \ref
   * mesh_functions are not well defined on points, e.g. the gradient of the
   * solution of a BVP is not well defined on the points of the mesh.
   * In this case WriteCellData(const std::string &name, const <!--
   * --> MESH_FUNCTION &mesh_function) may be more appropriate.

   *
   * ### Example usage
   * @snippet vtk_writer.cc mfPointUsage
   */
  template <class MESH_FUNCTION,
            class = std::enable_if_t<uscalfe::isMeshFunction<MESH_FUNCTION>>>
  void WritePointData(const std::string& name,
                      const MESH_FUNCTION& mesh_function);

  /**
   * @brief Add a new `unsigned char` attribute dataset that attaches data to
   * the `cells` of the mesh (i.e. to entities with codim = "codim that was
   * specified in constructor of VtkWriter")
   * @param name The name of the attribute set.
   * @param mds The mesh dataset that that attaches the data to the cells
   *            the mesh.
   * @param undefined_value The value that should be written for a cell to
   * which `mds` does not attach data (i.e. if `mds.DefinedOn() == false`)
   */
  void WriteCellData(const std::string& name,
                     const mesh::utils::MeshDataSet<unsigned char>& mds,
                     unsigned char undefined_value = 0);

  /**
   * @brief Add a new `char` attribute dataset that attaches data to
   * the `cells` of the mesh (i.e. to entities with codim = "codim that was
   * specified in constructor of VtkWriter")
   * @param name The name of the attribute set.
   * @param mds The mesh dataset that that attaches the data to the cells
   *            the mesh.
   * @param undefined_value The value that should be written for a cell to
   * which `mds` does not attach data (i.e. if `mds.DefinedOn() == false`)
   */
  void WriteCellData(const std::string& name,
                     const mesh::utils::MeshDataSet<char>& mds,
                     char undefined_value = 0);

  /**
   * @brief Add a new `unsigned int` attribute dataset that attaches data to
   * the `cells` of the mesh (i.e. to entities with codim = "codim that was
   * specified in constructor of VtkWriter")
   * @param name The name of the attribute set.
   * @param mds The mesh dataset that that attaches the data to the cells
   *            the mesh.
   * @param undefined_value The value that should be written for a cell to
   * which `mds` does not attach data (i.e. if `mds.DefinedOn() == false`)
   */
  void WriteCellData(const std::string& name,
                     const mesh::utils::MeshDataSet<unsigned int>& mds,
                     unsigned int undefined_value = 0);

  /**
   * @brief Add a new `int` attribute dataset that attaches data to
   * the `cells` of the mesh (i.e. to entities with codim = "codim that was
   * specified in constructor of VtkWriter")
   * @param name The name of the attribute set.
   * @param mds The mesh dataset that that attaches the data to the cells
   *            the mesh.
   * @param undefined_value The value that should be written for a cell to
   * which `mds` does not attach data (i.e. if `mds.DefinedOn() == false`)
   */
  void WriteCellData(const std::string& name,
                     const mesh::utils::MeshDataSet<int>& mds,
                     int undefined_value = 0);

  /**
   * @brief Add a new `float` attribute dataset that attaches data to
   * the `cells` of the mesh (i.e. to entities with codim = "codim that was
   * specified in constructor of VtkWriter")
   * @param name The name of the attribute set.
   * @param mds The mesh dataset that that attaches the data to the cells
   *            the mesh.
   * @param undefined_value The value that should be written for a cell to
   * which `mds` does not attach data (i.e. if `mds.DefinedOn() == false`)
   */
  void WriteCellData(const std::string& name,
                     const mesh::utils::MeshDataSet<float>& mds,
                     float undefined_value = 0);

  /**
   * @brief Add a new `double` attribute dataset that attaches data to
   * the `cells` of the mesh (i.e. to entities with codim = "codim that was
   * specified in constructor of VtkWriter")
   * @param name The name of the attribute set.
   * @param mds The mesh dataset that that attaches the data to the cells
   *            the mesh.
   * @param undefined_value The value that should be written for a cell to
   * which `mds` does not attach data (i.e. if `mds.DefinedOn() == false`)
   */
  void WriteCellData(const std::string& name,
                     const mesh::utils::MeshDataSet<double>& mds,
                     double undefined_value = 0);

  /**
   * @brief Add a new vector attribute dataset that attaches vectors to
   * the cells of the mesh.
   * @param name The name of the attribute set.
   * @param mds The mesh dataset that that attaches the data to the cells of
   *            the mesh.
   * @param undefined_value The value that should be written for a cell to
   * which `mds` does not attach data (i.e. if `mds.DefinedOn() == false`)
   */
  void WriteCellData(const std::string& name,
                     const mesh::utils::MeshDataSet<Eigen::Vector2d>& mds,
                     const Eigen::Vector2d& undefined_value = {0, 0});

  /**
   * @brief Add a new vector attribute dataset that attaches vectors to
   * the cells of the mesh.
   * @param name The name of the attribute set.
   * @param mds The mesh dataset that that attaches the data to the cells of
   *            the mesh.
   * @param undefined_value The value that should be written for a cell to
   * which `mds` does not attach data (i.e. if `mds.DefinedOn() == false`)
   */
  void WriteCellData(const std::string& name,
                     const mesh::utils::MeshDataSet<Eigen::Vector2f>& mds,
                     const Eigen::Vector2f& undefined_value = {0, 0});

  /**
   * @brief Add a new vector attribute dataset that attaches vectors to
   * the cells of the mesh.
   * @param name The name of the attribute set.
   * @param mds The mesh dataset that that attaches the data to the cells of
   *            the mesh.
   * @param undefined_value The value that should be written for a cell to
   * which `mds` does not attach data (i.e. if `mds.DefinedOn() == false`)
   */
  void WriteCellData(const std::string& name,
                     const mesh::utils::MeshDataSet<Eigen::Vector3d>& mds,
                     const Eigen::Vector3d& undefined_value = {0, 0, 0});

  /**
   * @brief Add a new vector attribute dataset that attaches vectors to
   * the cells of the mesh.
   * @param name The name of the attribute set.
   * @param mds The mesh dataset that that attaches the data to the cells of
   *            the mesh.
   * @param undefined_value The value that should be written for a cell to
   * which `mds` does not attach data (i.e. if `mds.DefinedOn() == false`)
   */
  void WriteCellData(const std::string& name,
                     const mesh::utils::MeshDataSet<Eigen::Vector3f>& mds,
                     const Eigen::Vector3f& undefined_value = {0, 0, 0});

  /**
   * @brief Add a new vector attribute dataset that attaches vectors to
   * the cells of the mesh.
   * @param name The name of the attribute set.
   * @param mds The mesh dataset that that attaches the data to the cells of
   *            the mesh.
   * @param undefined_value The value that should be written for a cell to
   * which `mds` does not attach data (i.e. if `mds.DefinedOn() == false`)
   * @note This version accepts in principle arbitrary size double Vectors,
   * however only the first three components are visualized!
   */
  void WriteCellData(
      const std::string& name,
      const mesh::utils::MeshDataSet<Eigen::VectorXd>& mds,
      const Eigen::VectorXd& undefined_value = Eigen::Vector3d(0, 0, 0));

  /**
   * @brief Add a new vector attribute dataset that attaches vectors to
   * the cells of the mesh.
   * @param name The name of the attribute set.
   * @param mds The mesh dataset that that attaches the data to the cells of
   *            the mesh.
   * @param undefined_value The value that should be written for a cell to
   * which `mds` does not attach data (i.e. if `mds.DefinedOn() == false`)
   * @note This version accepts in principle arbitrary size float Vectors,
   * however only the first three components are visualized!
   */
  void WriteCellData(
      const std::string& name,
      const mesh::utils::MeshDataSet<Eigen::VectorXf>& mds,
      const Eigen::VectorXf& undefined_value = Eigen::Vector3f(0, 0, 0));

  /**
   * @brief Sample a \ref mesh_function "MeshFunction" at the barycenter of the
   * cell and visualize it as cell data in the Vtk File.
   * @tparam MESH_FUNCTION An object fulfilling the \ref mesh_function
   * concept. The \ref uscalfe::MeshFunctionReturnType should be one of
   * - unsigned char
   * - char
   * - unsigned int
   * - int
   * - float
   * - double
   * - Eigen::Vector2d
   * - Eigen::Vector2f
   * - Eigen::Vector3d
   * - Eigen::Vector3f
   * - Eigen::VectorXd
   * - Eigen::VectorXf
   * @param name The name of the dataset, shouldn't contain any spaces!
   * @param mesh_function The \ref mesh_function "Mesh Function" to be sampled.
   *
   * @note this function will evaluate the MeshFunction on barycenters
   * of the cells of the mesh, i.e. at codim=0 entities. Some \ref
   * mesh_functions are not well defined on cells (e.g. boundary conditions)
   * and cannot be visualized with this function. (An assert will fail)
   *
   * ### Example usage
   * @snippet vtk_writer.cc mfCellUsage
   */
  template <class MESH_FUNCTION,
            class = std::enable_if_t<uscalfe::isMeshFunction<MESH_FUNCTION>>>
  void WriteCellData(const std::string& name,
                     const MESH_FUNCTION& mesh_function);

  /**
   * @brief Write global data into the vtk file that is not related to the mesh
   *        at all.
   * @param name Name of the global dataset.
   * @param data The data series that should be written.
   *
   * This function writes so-called "Field Data" into the vtk file which can
   * be read and visualized by paraview. It is e.g. a convenient way to export
   * the simulation time or the global mesh size.
   */
  void WriteGlobalData(const std::string& name, std::vector<int> data);

  /**
   * @brief Write global data into the vtk file that is not related to the mesh
   *        at all.
   * @param name Name of the global dataset.
   * @param data The data series that should be written.
   *
   * This function writes so-called "Field Data" into the vtk file which can
   * be read and visualized by paraview. It is e.g. a convenient way to export
   * the simulation time or the global mesh size.
   */
  void WriteGlobalData(const std::string& name, std::vector<float> data);

  /**
   * @brief Write global data into the vtk file that is not related to the mesh
   *        at all.
   * @param name Name of the global dataset.
   * @param data The data series that should be written.
   *
   * This function writes so-called "Field Data" into the vtk file which can
   * be read and visualized by paraview. It is e.g. a convenient way to export
   * the simulation time or the global mesh size.
   */
  void WriteGlobalData(const std::string& name, std::vector<double> data);

  /**
   * @brief Destructor, writes everything into the file and closes it.
   */
  ~VtkWriter() { WriteToFile(vtk_file_, filename_); }

 private:
  std::shared_ptr<const mesh::Mesh> mesh_;
  VtkFile vtk_file_;
  std::string filename_;
  dim_t codim_;
  unsigned char order_;

  // entry [i] stores the reference coordinates for the auxilliary nodes for
  // RefEl.Id() == i
  std::array<Eigen::MatrixXd, 5> aux_nodes_;

  // entry[cd] contains the index of the first auxiliary node for the codim=cd
  // entity
  std::array<mesh::utils::CodimMeshDataSet<unsigned int>, 2> aux_node_offset_;

  template <class T>
  void WriteScalarPointData(const std::string& name,
                            const mesh::utils::MeshDataSet<T>& mds,
                            T undefined_value);

  template <int ROWS, class T>
  void WriteVectorPointData(
      const std::string& name,
      const mesh::utils::MeshDataSet<Eigen::Matrix<T, ROWS, 1>>& mds,
      const Eigen::Matrix<T, ROWS, 1>& undefined_value);

  template <class T>
  void WriteScalarCellData(const std::string& name,
                           const mesh::utils::MeshDataSet<T>& mds,
                           T undefined_value);

  template <int ROWS, class T>
  void WriteVectorCellData(
      const std::string& name,
      const mesh::utils::MeshDataSet<Eigen::Matrix<T, ROWS, 1>>& mds,
      const Eigen::Matrix<T, ROWS, 1>& undefined_value);

  template <class T>
  void WriteFieldData(const std::string& name, std::vector<T> data);

  template <class DATA>
  void CheckAttributeSetName(const DATA& data, const std::string& name);

  template <int ROWS, class T>
  static void PadWithZeros(Eigen::Matrix<T, 3, 1>& out,
                           const Eigen::Matrix<T, ROWS, 1>& in);
};

template <class MESH_FUNCTION, class>
void VtkWriter::WritePointData(const std::string& name,
                               const MESH_FUNCTION& mesh_function) {
  if (order_ == 1) {
    Eigen::Matrix<double, 0, 1> origin{};
    WritePointData(
        name, *mesh::utils::make_LambdaMeshDataSet(
                  [&](const auto& e) { return mesh_function(e, origin)[0]; }));
  } else {
    // for higher orders, we have to evaluate the mesh function at all points
    // and write this into the vtk file:
    CheckAttributeSetName(vtk_file_.point_data, name);
    using T = uscalfe::MeshFunctionReturnType<MESH_FUNCTION>;
    auto dim_mesh = mesh_->DimMesh();

    Eigen::Matrix<double, 0, 1> origin;

    if constexpr (std::is_same_v<T, unsigned char> || std::is_same_v<T, char> ||
                  std::is_same_v<T, unsigned> || std::is_same_v<T, int> ||
                  std::is_same_v<T, float> || std::is_same_v<T, double>) {
      // MeshFunction is scalar valued:
      VtkFile::ScalarData<T> data{};
      data.data.resize(vtk_file_.unstructured_grid.points.size());
      data.name = name;

      // evaluate at nodes:
      for (auto& n : mesh_->Entities(dim_mesh)) {
        data.data[mesh_->Index(n)] = mesh_function(n, origin)[0];
      }

      for (char codim = dim_mesh - 1; codim >= static_cast<char>(codim_);
           --codim) {
        for (auto& e : mesh_->Entities(codim)) {
          auto ref_el = e.RefEl();
          if (order_ < 3 && ref_el == base::RefEl::kTria()) {
            continue;
          }
          auto values = mesh_function(e, aux_nodes_[ref_el.Id()]);
          auto offset = aux_node_offset_[codim](e);
          for (int i = 0; i < values.size(); ++i) {
            data.data[offset + i] = values[i];
          }
        }
      }
      vtk_file_.point_data.push_back(std::move(data));
    } else if constexpr (base::is_eigen_matrix<T>) {
      static_assert(T::ColsAtCompileTime == 1,
                    "The MeshFunction must return row-vectors");
      static_assert(
          T::RowsAtCompileTime == Eigen::Dynamic ||
              (T::RowsAtCompileTime > 1 && T::RowsAtCompileTime < 4),
          "The Row vectors returned by the MeshFunction must be either dynamic "
          "or must have 2 or 3 rows (at compile time).");
      using Scalar = typename T::Scalar;
      static_assert(
          std::is_same_v<double, Scalar> || std::is_same_v<float, Scalar>,
          "The RowVectors returned by the MeshFunction must be either double "
          "or float valued.");
      VtkFile::VectorData<Scalar> data{};
      data.data.resize(vtk_file_.unstructured_grid.points.size());
      data.name = name;

      // evaluate at nodes:
      for (auto& n : mesh_->Entities(dim_mesh)) {
        PadWithZeros<T::RowsAtCompileTime, Scalar>(data.data[mesh_->Index(n)],
                                                   mesh_function(n, origin)[0]);
      }

      for (char codim = dim_mesh - 1; codim >= static_cast<char>(codim_);
           --codim) {
        for (auto& e : mesh_->Entities(codim)) {
          auto ref_el = e.RefEl();
          if (order_ < 3 && ref_el == base::RefEl::kTria()) {
            continue;
          }
          auto values = mesh_function(e, aux_nodes_[ref_el.Id()]);
          auto offset = aux_node_offset_[codim](e);
          for (int i = 0; i < values.size(); ++i) {
            PadWithZeros<T::RowsAtCompileTime, Scalar>(data.data[offset + i],
                                                       values[i]);
          }
        }
      }
      vtk_file_.point_data.push_back(std::move(data));
    } else {
      LF_VERIFY_MSG(false,
                    "MeshFunction values must be one of: unsigned char, char, "
                    "unsigned, int, float, double, Eigen::Vector<double, ...> "
                    "or Eigen::Vector<float, ...>");
    }
  }
}

template <class MESH_FUNCTION, class>
void VtkWriter::WriteCellData(const std::string& name,
                              const MESH_FUNCTION& mesh_function) {
  // maps from RefEl::Id() -> barycenter of the reference element
  std::vector<Eigen::VectorXd> barycenters(5);
  barycenters[base::RefEl::kPoint().Id()] = Eigen::Matrix<double, 0, 1>();
  barycenters[base::RefEl::kSegment().Id()] =
      base::RefEl::kSegment().NodeCoords().rowwise().sum() / 2.;
  barycenters[base::RefEl::kTria().Id()] =
      base::RefEl::kTria().NodeCoords().rowwise().sum() / 3.;
  barycenters[base::RefEl::kQuad().Id()] =
      base::RefEl::kTria().NodeCoords().rowwise().sum() / 4.;

  WriteCellData(name, *mesh::utils::make_LambdaMeshDataSet([&](const auto& e) {
                  return mesh_function(e, barycenters[e.RefEl().Id()])[0];
                }));
}

template <class DATA>
void VtkWriter::CheckAttributeSetName(const DATA& data,
                                      const std::string& name) {
  if (std::find_if(data.begin(), data.end(), [&](auto& d) {
        return boost::apply_visitor([&](auto&& d2) { return d2.name; }, d) ==
               name;
      }) != data.end()) {
    throw base::LfException(
        "There is already another Point/Cell Attribute Set with the "
        "name " +
        name);
  }
  if (name.find(' ') != std::string::npos) {
    throw base::LfException(
        "The name of the attribute set cannot contain spaces!");
  }
}

template <int ROWS, class T>
void VtkWriter::PadWithZeros(Eigen::Matrix<T, 3, 1>& out,
                             const Eigen::Matrix<T, ROWS, 1>& in) {
  if constexpr (ROWS == 2) {  // NOLINT
    out.template block<2, 1>(0, 0) = in;
    out(2) = T(0);
  } else if constexpr (ROWS == 3) {  // NOLINT
    out = in;
  } else if constexpr (ROWS == Eigen::Dynamic) {  // NOLINT
    if (in.rows() == 2) {
      out(2) = T(0);
      out.topRows(in.rows()) = in;
    } else {
      out = in;
    }
  }
}

}  // namespace lf::io

#endif  // __3e48c7b32a034cb3be3dbca884ff4f6c
