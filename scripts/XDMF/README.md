# XDMF writer usage and dependency setup

## Dependency overview

The following libraries are needed for the XDMF writer / checkpointing system to function properly:

- HDF5 v1.10.4+
  - MPI parallelism (`--enable-parallel`)
  - C++ HL bindings (`--enable-hl`), required by ASAGI
  - Fortran bindings (`--enable-fortran`)
  - Support for multithreaded usage (`--enable-threadsafe`, might require `--enable-unsupported`)

- FoX
  - Support for dynamic linking (`-DBUILD_SHARED_LIBS=1`), the author provides a GitHub fork where this option is implemented
  - Suppress non-fatal, but annoying warnings (`-DFoX_SUPPRESS_WARNINGS=1`), same fork as above
  - XML DOM support (`-DFoX_ENABLE_DOM=1`)
  - Disable unused modules (`-DFoX_ENABLE_EXAMPLES=0 -DFoX_ENABLE_WCML=0 -DFoX_ENABLE_WKML=0 -DFoX_ENABLE_WXML=0`)

When the ASAGI server is used, ASAGI has to be build and linked against the same HDF5 libraries as Samoa. This is to prevent library conflicts or data corruption. Thus, the dependencies of ASAGI have to be compiled and linked against the correct HDF5 libraries.

- NetCDF C library
  - MPI parallelism (`--enable-parallel4`)
  - Disable networking, this removes the libcurl dependency (`--disable-dap`)

- NetCDF C++ library

### Binaries

In the end, you should have the following libraries in your link path:
- `libhdf5.so`
- `libhdf5_fortran.so`
- `libFoX_dom.so`
- `libFoX_sax.so`
- `libFoX_wxml.so`
- `libFoX_common.so`
- `libFoX_fsys.so`
- `libFoX_utils.so`

And for ASAGI in addition:
- `libasagi.so`
- `libhdf5hl_fortran.so`
- `libhdf5_hl.so`
- `libnetcdf.so`

Note that additional static libraries (`.a`) may be present, but these are not used at the time writing.
For a list of files that should be in your include path, refer to appendix A.

## Custom built vs. pre-installed packages

Some Linux distributions do not yet package the recent required versions of HDF5 and FoX. In such cases the compilation from source is recommended.

Using the the pre-installed modules on CoolMUC and SuperMUC-NG is possible. The recommended modules are:
- hdf5/1.10.2-intel-impi-threadsafe

For ASAGI:
- netcdf/4.6.1-intel-impi-hdf5v1.10-parallel
- netcdf-cxx4/4.3.0-intel-impi-hdf5v1.10
- netcdf-fortran/4.4.4-intel-impi-hdf5v1.10

## Library installation scripts

To ensure the correct compilation und installation of all needed dependencies, a set of `bash` shell installation scripts are provided in this directory. The scripts are able to download and install the `hdf5`, `fox`, `netcdf_c`, `netcdf_cpp` and `asagi` dependencies.

**Important!** The default and recommended behavior is to install all libraries into a directory outside of the default system search path. This is to prevent library conflicts with the existing system. However, if you already have a version of HDF5 and/or FoX, ASAGI or NetCDF in your search path (either by installing from the distribution package sources or by loading a module), this might cause trouble. Check your link path and order carefully.

The scripts require `git`, `autotools` and `cmake` in addition to the GNU compiler or intel compiler toolchain.

The scripts can be run in 4 modes of configuration:
 - GNU compiler without MPI `gnu nompi`
 - GNU compiler with MPI `gnu mpi`
 - Intel compiler without MPI `intel nompi`
 - Intel compiler with MPI `intel mpi`

The following environment variables may be configured:
 - `GIT_PROTOCOL`. May be either `https` (default) or `ssh`. Use this if you have to tunnel your outgoing internet connections (e.g. on SuperMUC-NG).
 - `ASAGI_NUMA`. Controls ASAGI NUMA support. Default: `-DNONUMA=1` (disabled NUMA). Set to `-DNONUMA=0` to enable NUMA support.
 - `HDF5_URL`. May be used to override the HDF5 git repository URL. Use this if you have no SSH access to the original HDF5 BitBucket repository.

### Scripts

- `./install_all_hpc.sh [library directory path]`: Install all required libraries into the specified directory, using the intel compiler with MPI support. If no directory is specified, the default location `~/local` will be used. It is advised to provide an absolute path, however testing whether this is strictly necessary is yet inconclusive. Use this script to setup your HPC environment.
- `./install_all.sh <mpi|nompi> <intel|gnu> [library directory path]`: This script is called by the `./install_all_hpc.sh` script with the parameters `mpi intel`. It can be used to specify the alternative GNU compiler or to specify the MPI support e.g. for testing.
- `./install_lib.sh <name> <mpi|nompi> <intel|gnu> [library directory path]`: This script is repeatedly called by the `./install_all.sh` script. It downloads and installs a specific library in the given configuration.
- `./asagi.sh`, `./fox.sh`, `./hdf5.sh`, `./netcdf_c.sh`, `./netcdf_cxx.sh`: These scripts are not meant to be called directly. Instead, use on of the wrapper scripts above. These scripts contain download and build instructions for the specific libraries. You may edit these in order to adjust the download sources or compilation flags.

After a successful compilation, the specified directory will contain the required libraries in `<directory>/<compiler>/<serial|parallel>/lib` and the header files in `<directory>/<compiler>/<serial|parallel>/include`. Feed this to your linker. For example, if you were to install the libraries to `/opt/samoa_xdmf_libs` using the intel compiler with MPI support, then the libraries would be located at `/opt/samoa_xdmf_libs/intel/parallel/lib`.

### Examples

- Local testing, using GCC and no MPI: `./install_all.sh nompi gnu /opt/samoa_xdmf_libs`
- Setup on SuperMUC-NG: `GIT_PROTOCOL=ssh ASAGI_NUMA='-DNONUMA=0' HDF5_URL='ssh://gitlab.lrz.de/ChristophHonal/hdf5.git' ./install_all_hpc.sh /dss/dsshome1/0D/ga63yos3/samoa_xdmf_libs`
  - This command requires a working SSH tunnel to GitHub and LRZ GitLab.
  - This command overrides the HDF5 repository URL, because the author has no SSH access to the HDF5 bitbucket, and SuperMUC-NG requires all outgoing connections to be via a SSH tunnel. The URL given is accessible to anyone logged in into LRZ GitLab.
  - This command installs the libraries to `/dss/dsshome1/0D/ga63yos3/samoa_xdmf_libs`, which is in the home directory of the author. Please adjust this to your own.

## Samoa configuration

Provided you are using a branch with proper XDMF support, the Samoa SCons configuration takes the following arguments, for example if the libraries were installed to `/opt/samoa_xdmf_libs`, using the GNU compiler with MPI support:

- `xdmf='true'`
- `xdmf_fox_dir='/opt/samoa_xdmf_libs/gnu/parallel'`
- `xdmf_hdf5_dir='/opt/samoa_xdmf_libs/gnu/parallel'`

When using ASAGI (`asagi='true'`), you must also configure the ASAGI and NetCDF linker paths:

- `asagi_dir='/opt/samoa_xdmf_libs/gnu/parallel'`
- `netcdf_dir='/opt/samoa_xdmf_libs/gnu/parallel'`


## Appendix A. List of files in include path

Note: This serves only as a rough reference, not as a must-have list of files. Your configuarion or installation may differ.

- `asagi.f90`
- `H5Apublic.h`
- `h5f.mod`
- `H5overflow.h`
- `hdf5.mod`
- `m_sax_reader.mod`
- `asagi.h`
- `H5Cpublic.h`
- `h5fortkit.mod`
- `H5PLextern.h`
- `m_common_attrs.mod`
- `m_sax_tokenizer.mod`
- `fox_common.mod`
- `h5d.mod`
- `h5fortran_types.mod`
- `H5PLpublic.h`
- `m_common_buffer.mod`
- `m_sax_types.mod`
- `fox_dom.mod`
- `H5DOpublic.h`
- `H5Fpublic.h`
- `h5p.mod`
- `m_common_charset.mod`
- `m_sax_xml_source.mod`
- `fox_m_fsys_abort_flush.mod`
- `H5Dpublic.h`
- `h5_gen.mod`
- `H5Ppublic.h`
- `m_common_content_model.mod`
- `m_wxml_core.mod`
- `fox_m_fsys_array_str.mod`
- `h5ds.mod`
- `h5global.mod`
- `H5PTpublic.h`
- `m_common_element.mod`
- `m_wxml_escape.mod`
- `fox_m_fsys_count_parse_input.mod`
- `H5DSpublic.h`
- `h5g.mod`
- `H5pubconf.h`
- `m_common_elstack.mod`
- `m_wxml_overloads.mod`
- `fox_m_fsys_format.mod`
- `h5e.mod`
- `H5Gpublic.h`
- `H5public.h`
- `m_common_entities.mod`
- `netcdf_aux.h`
- `fox_m_fsys_parse_input.mod`
- `H5Epubgen.h`
- `h5im.mod`
- `h5r.mod`
- `m_common_entity_expand.mod`
- `netcdf_dispatch.h`
- `fox_m_fsys_realtypes.mod`
- `H5Epublic.h`
- `h5i.mod`
- `H5Rpublic.h`
- `m_common_error.mod`
- `netcdf_filter.h`
- `fox_m_fsys_string_list.mod`
- `H5f90i_gen.h`
- `H5IMpublic.h`
- `h5s.mod`
- `m_common_io.mod`
- `netcdf.h`
- `fox_m_fsys_string.mod`
- `H5f90i.h`
- `H5Ipublic.h`
- `H5Spublic.h`
- `m_common_namecheck.mod`
- `netcdf_mem.h`
- `fox_m_fsys_varstr.mod`
- `H5FDcore.h`
- `H5LDpublic.h`
- `h5tb_const.mod`
- `m_common_namespaces.mod`
- `netcdf_meta.h`
- `fox_m_utils_mtprng.mod`
- `H5FDdirect.h`
- `h5lib.mod`
- `h5tb.mod`
- `m_common_notations.mod`
- `netcdf_par.h`
- `fox_m_utils_uri.mod`
- `H5FDfamily.h`
- `h5l.mod`
- `H5TBpublic.h`
- `m_common_struct.mod`
- `fox_m_utils_uuid.mod`
- `H5FDlog.h`
- `H5Lpublic.h`
- `h5t.mod`
- `m_dom_dom.mod`
- `fox_sax.mod`
- `H5FDmpi.h`
- `h5lt_const.mod`
- `H5Tpublic.h`
- `m_dom_error.mod`
- `fox_utils.mod`
- `H5FDmpio.h`
- `h5lt.mod`
- `H5version.h`
- `m_dom_extras.mod`
- `fox_wxml.mod`
- `H5FDmulti.h`
- `H5LTpublic.h`
- `h5z.mod`
- `m_dom_parse.mod`
- `H5ACpublic.h`
- `H5FDpublic.h`
- `H5MMpublic.h`
- `H5Zpublic.h`
- `m_dom_utils.mod`
- `h5a.mod`
- `H5FDsec2.h`
- `h5o.mod`
- `hdf5.h`
- `m_sax_operate.mod`
- `H5api_adpt.h`
- `H5FDstdio.h`
- `H5Opublic.h`
- `hdf5_hl.h`
- `m_sax_parser.mod`