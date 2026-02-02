block_cipher = None
import site
import sys
sys.setrecursionlimit(sys.getrecursionlimit() * 5)
p = site.getsitepackages()[-1] + '/'



from PyInstaller.utils.hooks import collect_data_files, collect_dynamic_libs

hiddenimports66 = [
    "pyarrow._parquet",
    "pyarrow.formatting",
    "pyarrow.lib",
    "pyarrow.vendored",
    "pyarrow.vendored.version",
    "pyarrow.compat",
    "pyarrow.compat.*"
]

datas66 = collect_data_files('pyarrow')

binaries = collect_dynamic_libs('pyarrow')

from PyInstaller.utils.hooks import collect_data_files

datas22 = collect_data_files('bokeh.core') + \
        collect_data_files('bokeh.server') + \
        collect_data_files('bokeh.command.subcommands', include_py_files=True) + \
        collect_data_files('bokeh')


hiddenimports = [
  'numpy',
  'dendropy',
  'scipy.integrate.lsoda',
  'sklearn.utils._cython_blas',
  'sklearn.neighbors.typedefs',
  'sklearn.neighbors.quad_tree',
  'sklearn.tree',
  'sklearn.tree._utils',
  'pandas._libs.tslibs.timedeltas',
  'scipy.special._ufuncs_cxx',
  'scipy.linalg.cython_blas',
  'scipy.linalg.cython_lapack',
  'scipy.integrate',
  'scipy.integrate.quadrature',
  'scipy.integrate.odepack',
  'scipy.integrate._odepack',
  'scipy.integrate.quadpack',
  'scipy.integrate._quadpack',
  'scipy.integrate._ode',
  'scipy.integrate.vode',
  'scipy.integrate._dop',
  'scipy.integrate.lsoda',
  'scipy._lib.messagestream' ,
] + hiddenimports66


a = Analysis(['gui.py'],
             binaries= binaries,
             datas=([('icon','icon'),
                     ('src','src'),

                   ('ccs.ui', '.'),
                   ('file_parser.py', '.'),
                   ('ccs_database.db', '.'),
                   ('popup.py', '.'),
                   ('psi-ms-2.5.0.obo.gz', '.'),
                   ('photo.qrc', '.'),
                   ('photo_rc.py', '.'),
                ] + datas22 + datas66),
             hiddenimports=hiddenimports,
             hookspath=[],
             runtime_hooks=[],
             excludes=[],
             win_no_prefer_redirects=False,
             win_private_assemblies=False,
             cipher=block_cipher,
             noarchive=False)
pyz = PYZ(a.pure, a.zipped_data,
             cipher=block_cipher)
exe = EXE(pyz,
         a.scripts,
         a.binaries,
         a.zipfiles,
         a.datas,
         [],
         name='CCS_find',
         debug=False,
         bootloader_ignore_signals=False,
         strip=False,
         upx=True,
         upx_exclude=[],
         runtime_tmpdir=None,
         console=False,
         icon='icon/ccs_logo.ico' )


coll = COLLECT(exe,
              a.binaries,
              a.zipfiles,
              a.datas,
              strip=False,
              upx=True,
              upx_exclude=[],
              name='CCS_find')
