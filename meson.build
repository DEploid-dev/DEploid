project('DEploid', ['cpp','c'],
        version: files('VERSION'),
        default_options : [
          'warning_level=3',
          'cpp_std=c++17'
        ],
        meson_version: '>=0.57')

cpp = meson.get_compiler('cpp')

# The problem here is that we don't know when to update these vars.
deploid_version = run_command( meson.current_source_dir() / 'deploid_version.sh' ).stdout().strip()
lassoversion = run_command( meson.current_source_dir() / 'lassoversion.sh' ).stdout().strip()
# If we are ALWAYS updating this one, then we will ALWAYS recompile!
# That would maybe be OK if we only do this in a specific small file.
compiledate = run_command( meson.current_source_dir() / 'compiledate.sh' ).stdout().strip()


add_project_arguments(['-DVERSION="'+meson.project_version() + '"',
                       '-DDEPLOID_VERSION="'+deploid_version+'"',
                       '-DLASSOVERSION="'+lassoversion+'"',
                       '-DCOMPILEDATE="'+compiledate+'"'],
                      language : 'cpp')

common_sources = ['src/random/fastfunc.cpp',
		  'src/random/random_generator.cpp',
		  'src/random/mersenne_twister.cpp',
		  'src/chooseK.cpp',
		  'src/ibd.cpp',
		  'src/mcmc.cpp',
		  'src/panel.cpp',
		  'src/utility.cpp',
		  'src/dEploidIO.cpp',
		  'src/dEploidIO-operation.cpp',
		  'src/dEploidIO-workflow.cpp',
		  'src/updateHap.cpp',
		  'src/txtReader.cpp',
		  'src/vcfReader.cpp',
		  'src/variantIndex.cpp',
		  'src/gzstream/gzstream.cpp',
		  'src/export/dEploidIOExport.cpp',
		  'src/export/dEploidIOExportPosteriorProb.cpp',
		  'src/export/writeMcmcRelated.cpp',
		  'src/lasso/src/lasso.cpp',
                  'src/lasso/src/dEploidLasso.cpp',
                  'src/debug/mcmcDebug.cpp',
                  'src/debug/vcfReaderDebug.cpp',
                  'src/lasso/src/lassoDBG.cpp']

src = include_directories('src')
gzstream_dir = include_directories('src/gzstream')
zlib = dependency('zlib')

libcommon = static_library('common', common_sources,
                           include_directories: [src, gzstream_dir],
                           dependencies: zlib)


executable('dEploid',
           ['src/dEploid.cpp'] + common_sources,
           include_directories: [src,gzstream_dir],
           dependencies: zlib,
           link_with: libcommon,
           install: true)

cppunit = dependency('cppunit')
libdl = cpp.find_library('dl')

unit_test_sources = {
                      'basic':   ['tests/unittest/test_ibd.cpp',
		                  'tests/unittest/test_updateSingleHap.cpp',
		                  'tests/unittest/test_utilities.cpp',
		                  'tests/unittest/test_panel.cpp',
		                  'tests/unittest/test_mcmc.cpp',
		                  'tests/unittest/test_updatePairHap.cpp',
		                  'tests/unittest/test_workflow.cpp'],
                      'io':      ['tests/unittest/test_dEploidIO.cpp',
			          'tests/unittest/test_txtReader.cpp',
			          'tests/unittest/test_vcfReader.cpp'],
                      'current': [ 'tests/unittest/test_workflow.cpp']
                    }

foreach test_name : unit_test_sources.keys()

  name = test_name+'-unit_tests'
  exe = executable( name,
                    ['tests/unittest/test_runner.cpp'] + unit_test_sources[test_name],
                    cpp_args: ['-DNDEBUG','-DUNITTEST','-Wno-write-strings','--coverage'],
                    link_args: ['--coverage'],
                    include_directories: [src, gzstream_dir],
                    dependencies: [zlib,cppunit,libdl],
                    link_with: libcommon)

  test(name + ' unit tests',
       exe,
       timeout: 200,
     workdir: meson.current_source_dir())

endforeach
     
