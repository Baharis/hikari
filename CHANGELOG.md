# CHANGELOG


## v0.3.2 (2025-01-22)

### :bug:

- :bug: Fix heuristic `Operation.origin` code that gave wrong results in cubic system
  ([`8ad71d1`](https://github.com/Baharis/hikari/commit/8ad71d18f12af13f25dcb80516fac85012a09961))

- :bug: Fix some trigonal systems being mistaken for cubic
  ([`51e50d1`](https://github.com/Baharis/hikari/commit/51e50d15dfb3d662e4f3ae5a8d7c96bc5e730d1d))

- :bug: Fix test which wrongly assumed that diagonal 3-fold screw axes should extinct -hhh, h-hh,
  and hh-h reflections
  ([`38bff51`](https://github.com/Baharis/hikari/commit/38bff519171f8eb9ab76bfc4c46a8285be097893))

- :bug: Fix test which wrongly assumed that diagonal 3-fold screw axes should extinct -hhh, h-hh,
  and hh-h reflections
  ([`6b98f10`](https://github.com/Baharis/hikari/commit/6b98f109f5bd6228158e014f851bcc836334fcfb))

- :bug: Make `Group.auto_generated_code` a single word
  ([`67c01c7`](https://github.com/Baharis/hikari/commit/67c01c75b01d49dff91e4566082db59ff19f85fe))

- :bug: Properly attach `Group.name` and `.number` when reading catalog from json
  ([`1952454`](https://github.com/Baharis/hikari/commit/1952454934b07729a89d606234444d630be6666d))

- :bug: Split `SymmOp` into `Operation` and `BoundedOperation` whose translation is fixed to [0,1)
  ([`5352e54`](https://github.com/Baharis/hikari/commit/5352e5448fea7ce420bdfd1a134130d2cbfa6b67))

- :bug: The direction of screw axes should align with their glide
  ([`ea7c587`](https://github.com/Baharis/hikari/commit/ea7c587081e05cac8c3061f2bfacdf39d069eeba))

### :rewind:

- :rewind: Revert ":test_tube: Add failing test for auto-generated space group names"
  ([`0f0cb4c`](https://github.com/Baharis/hikari/commit/0f0cb4ce17d3da438980a49f639752f5fdf4c963))

This reverts commit 50c96c21928f208b7834e54e25f2979885e08ad8.

### Other

- :bulb: Remove TODO comment associated with completed code update
  ([`80189b8`](https://github.com/Baharis/hikari/commit/80189b80579d2a27a926df1251c1f5af3b5c6265))

- :construction: Add `hikari.utility.Singleton` to be potentially used later
  ([`a051311`](https://github.com/Baharis/hikari/commit/a051311147cd8639538a1f3e3f4fba20f163b85e))

- :construction: Add rudimentary `hm_symbol` property that better aligns with crystallographic
  tables
  ([`8b850b0`](https://github.com/Baharis/hikari/commit/8b850b0df6eb833a312e7dcad40fb44c082c067f))

- :memo: Improve documentation, type hints in `hikari.dataframes.cif`
  ([`2018a6c`](https://github.com/Baharis/hikari/commit/2018a6c2a5d834e4432aa410cfaa34cdf8740862))

- :memo: Improve documentation, type hints in `hikari.dataframes.cif`
  ([`941c6a1`](https://github.com/Baharis/hikari/commit/941c6a1bc0439c45775e10219ac2412641b807c2))

- :memo: Improve type hints for `hikari.symmetry.Operation`s
  ([`ce974d0`](https://github.com/Baharis/hikari/commit/ce974d0aaec811de731979dfd6ea89438cf9151c))

- :memo: Streamline `operations.py` type hints to rely on `typing` only
  ([`3d003ab`](https://github.com/Baharis/hikari/commit/3d003ab8b9ed3821f1eee70690926f127fab9b43))

- :test_tube: Add failing test for auto-generated space group names
  ([`50c96c2`](https://github.com/Baharis/hikari/commit/50c96c21928f208b7834e54e25f2979885e08ad8))

- :white_check_mark: Make (ub)bounded a property, add/improve all `Operation` tests
  ([`ad33a95`](https://github.com/Baharis/hikari/commit/ad33a95a76308ebd5bea35a5f352a634d18490c7))

- Merge remote-tracking branch 'origin/development' into development
  ([`321bc56`](https://github.com/Baharis/hikari/commit/321bc56b647e7987824579a6fd79636e3080e364))

- 🔀 Merge pull request #43 from Baharis/development
  ([`c2f66e9`](https://github.com/Baharis/hikari/commit/c2f66e98e6c3bb9270d0d0f904c5cc0fdc2e9dbf))

Split `SymmOp` into `Operation` and `BoundedOperation`


## v0.3.1 (2024-12-10)

### :bug:

- :bug: Fix `GroupCatalog.values`, `.items`
  ([`81d700c`](https://github.com/Baharis/hikari/commit/81d700c9ca1ac453b42ac6c85d30fa00c2e1a8f1))

- :bug: Fix all space groups being "standard" due to wrong type in json
  ([`ebe9b7d`](https://github.com/Baharis/hikari/commit/ebe9b7db6b804bea2fadcac85381f6b1babca36b))

### :package:

- :package: Add `MANIFEST.in` to potentially help adding resources to wheel
  ([`2875cf7`](https://github.com/Baharis/hikari/commit/2875cf70dd366270738aea5aa947f197c1cad4cb))

### Other

- :coffin: Remove dead resources function `_load_bytes`
  ([`de35f29`](https://github.com/Baharis/hikari/commit/de35f29f92c3f6be0f4560cbc96ea488a7510c4c))

- :construction_worker: In continuous integration, run rapid tests only
  ([`49bc6a9`](https://github.com/Baharis/hikari/commit/49bc6a977a4c677a2347dc8f4c86229d0e359e3a))

- :green_heart: Fix continuous integration to properly publish `hikari` under name `hikari-toolkit`
  ([`3992db7`](https://github.com/Baharis/hikari/commit/3992db795fcb2d4c414db8914849eb55dff0781f))

- :white_check_mark: Add even more tests for `GroupCatalog`
  ([`b6b17d6`](https://github.com/Baharis/hikari/commit/b6b17d6a0a4bc7391afa2a557f3287d92e6018e4))

- :white_check_mark: Add further tests for `GroupCatalog`
  ([`1723078`](https://github.com/Baharis/hikari/commit/17230789b0134be3b4fc8db3d3edaa52ea322254))

- :white_check_mark: Skip the potency map tests on python 3.13.0 only
  ([`b8f178a`](https://github.com/Baharis/hikari/commit/b8f178a816ec5be28d47d7a8f75fc61dce1fba30))

- 🔀 Merge pull request #42 from Baharis/development
  ([`fe1ca24`](https://github.com/Baharis/hikari/commit/fe1ca24d0220ac8ea69ef2a8739bcef2e06e4281))

Fix `GroupCatalog`, fix and expand `GroupCatalog` testing


## v0.3.0 (2024-12-02)

### :pencil2:

- :pencil2: Fix typo in pyproject.toml path
  ([`4dde1ef`](https://github.com/Baharis/hikari/commit/4dde1ef49ea0c718cf50796d66e31d281c7e75f4))

### :sparkles:

- :sparkles: Increase overall quality of handling, storing, selecting `Groups`, Increase user
  experience quality and code safety, Bug fixes.
  ([`ec6aac6`](https://github.com/Baharis/hikari/commit/ec6aac69f516c5bf812588fa7ec1d8cdf773bae1))

### Other

- .* suffix can only be used with `==` or `!=` operators
  ([`181004d`](https://github.com/Baharis/hikari/commit/181004d1515a0540b73b63a4d5d8d35a3dbc4fbd))

- :beers: Add forgotten type fix
  ([`30f65d7`](https://github.com/Baharis/hikari/commit/30f65d7da4fc7c5eb1abf7f1aec137a37396ae81))

- :green_heart: Add hikari API tokens
  ([`847188d`](https://github.com/Baharis/hikari/commit/847188de7f0f1353f624cb7c50319fe33399c096))

- :green_heart: Can't install scipy on newest macs, maybe issue with 3.13; skip for now
  ([`3a613bc`](https://github.com/Baharis/hikari/commit/3a613bc8f76913f99d27b52a3c2460560167215f))

- :green_heart: Make macos install gfortran
  ([`9d5da57`](https://github.com/Baharis/hikari/commit/9d5da5726754e774d53a6ad1827d5d353c73c915))

- :green_heart: Make macos install gfortran 2
  ([`85d3d9e`](https://github.com/Baharis/hikari/commit/85d3d9e08895a66dcdd4c058a18e12cc1b483ac0))

- :green_heart: Make macos install gfortran 3
  ([`b349f92`](https://github.com/Baharis/hikari/commit/b349f92d469c8f6ec0f3e3ef8ac0a0c3facd5d5b))

- :green_heart: Make macos install gfortran 4
  ([`32175c0`](https://github.com/Baharis/hikari/commit/32175c0fbcb05530e23fa2d6d78e745fd944cc23))

- :memo: Align documentation style with that of picometer
  ([`f25ab46`](https://github.com/Baharis/hikari/commit/f25ab467fc25e4912963720e8f5f16418854f27e))

- :tada: Move to poetry package control, bump required python version to 3.9
  ([`38397d1`](https://github.com/Baharis/hikari/commit/38397d1621a6336da13b191972915c4ded2d2832))

- :truck: Rename `test` directory to `tests`, remove some mentions of `picometer`
  ([`279926f`](https://github.com/Baharis/hikari/commit/279926f4c873ebb9792d4b4327b11f1dad3cf0b2))

- `groupcatalogue` to make group dicts, `hikari.typing`, `a`->`e`
  ([`491b1cb`](https://github.com/Baharis/hikari/commit/491b1cbd9a2569bcbd3782a3ea57166a72312d1e))

- Add a custom `AmbiguousGroupAccessorWarning` for clarity
  ([`8c8ab2e`](https://github.com/Baharis/hikari/commit/8c8ab2ea9bf6f92db2370bc857569fe43604185a))

- Add JSON files to git
  ([`414d047`](https://github.com/Baharis/hikari/commit/414d047fe7fb273dc949fd31f11460af5fc15db5))

- Add P21/n to the list of known structures (but not to pickles)
  ([`fa4d928`](https://github.com/Baharis/hikari/commit/fa4d928c26a94c645d1a3f05fcfa516a956c7f67))

- Add P21/n to the list of known structures2 (but not to pickles)
  ([`b624353`](https://github.com/Baharis/hikari/commit/b624353a0ba24b421f470312805f2f2770257e8b))

- Add python 3.10 to test builds
  ([`6d87228`](https://github.com/Baharis/hikari/commit/6d87228999adae46b945bcb320f3e88a8ee27162))

- Add some type hints
  ([`67d00de`](https://github.com/Baharis/hikari/commit/67d00dea5f28a3831ee3b6b205ecd6b8925734d9))

- Add tests for hkl completeness scripts
  ([`1026827`](https://github.com/Baharis/hikari/commit/1026827b03869c026e7d0dac8648f90b053a3509))

- Add tests for HklFrame, including writing hkl and res
  ([`cd2d59f`](https://github.com/Baharis/hikari/commit/cd2d59f840cbcb1cd791061d66124b7b3189b1c4))

- Add windows and mac os to test environments
  ([`e528da4`](https://github.com/Baharis/hikari/commit/e528da4c23261a7e3e80955a822da7034eb35bbe))

- Americanize all instances of `Catalogue` to `Catalog`
  ([`3438f82`](https://github.com/Baharis/hikari/commit/3438f82893122959e9efc06fedd1cd7e34934ebc))

- Apparently np.warnings doesn't exist in some versions; just try
  ([`f897412`](https://github.com/Baharis/hikari/commit/f8974127d40136856eee4977b53a9be098bfa8bb))

- Apparently np.warnings doesn't exist in some versions; just try
  ([`62b62e6`](https://github.com/Baharis/hikari/commit/62b62e6189d22211c325d40ffdfb64cbc7e3a38e))

- Apparently one refl is not always read correctly in Python3.10
  ([`9159136`](https://github.com/Baharis/hikari/commit/9159136706cd8f2f8b259f71fcef4f4c6bda1e74))

- Automatically generate the `GroupCatalogue` based on provided `GroupCatalogueKey` data
  ([`c6096a1`](https://github.com/Baharis/hikari/commit/c6096a1dccb30481d547fe011819f8074194269e))

- Automatically get correct space group based on index, HM, HM-short*, or Hall name
  ([`b9ee7f8`](https://github.com/Baharis/hikari/commit/b9ee7f86fca42ffa8c5a7da6f33045a2d56983e2))

- Bugfix: correctly type hint `Match`, remove fixed `TODO`s
  ([`fff5255`](https://github.com/Baharis/hikari/commit/fff525521a09a882188d4e4792e98eb65204c827))

- Compose a `Hall_symbol_PG.dat` for point groups only
  ([`6009f6d`](https://github.com/Baharis/hikari/commit/6009f6d9242f9864ffc2c913f899e727d1d12fa4))

- Define all keys using clearer OOP syntax with registry
  ([`98e79ee`](https://github.com/Baharis/hikari/commit/98e79eeb78ce2e7e74858e99f1f6a51539cab1bf))

- Define new `HklKeyDict` to substitute old ambiguous `HklKeys`
  ([`c06d0e3`](https://github.com/Baharis/hikari/commit/c06d0e3cb7291b6ae4438119a767aa2571aa5032))

- Do not use `>=` when specified starred version of dependency
  ([`c03b44a`](https://github.com/Baharis/hikari/commit/c03b44a106e9ba7f04a50cf4ce1d8831af354c71))

- Don't suggest buggy pandas 1.5.0 in `requirements.txt`
  ([`9908054`](https://github.com/Baharis/hikari/commit/9908054e359be3d67b467731a3c7aa8d4241f065))

- Except ax.dist deprecation; any alternatives?
  ([`1ea5270`](https://github.com/Baharis/hikari/commit/1ea52709e56822ee42ce71150d37b4462f9b3ecd))

- Fix `mpl` deprecation issues, `uncertainties` not supported...
  ([`7e37464`](https://github.com/Baharis/hikari/commit/7e374644defc3dab12f5cd4ccbdcee29081c5096))

- Fix code test environment to ubuntu-20.04
  ([`7764931`](https://github.com/Baharis/hikari/commit/7764931a48f5dfb78a1d2e94b0a7dd9dc1cf7785))

- Fix error in shelx_4 hkl format docstring
  ([`6c7ae78`](https://github.com/Baharis/hikari/commit/6c7ae78c70eb4b2d2cd46b0a5848d5a2cfac5cbc))

- Fix implementation for point groups
  ([`4fed404`](https://github.com/Baharis/hikari/commit/4fed4047023abbbe8efc9291c7802cbd6ef56c0a))

- Fix issues with new `HklKeyDict` and `HklKeyRegistrar`
  ([`fea9177`](https://github.com/Baharis/hikari/commit/fea917703da73af1da885af867578a4438e7f07c))

- Fix platform-specific `np.str` to `np.str_`
  ([`7400d26`](https://github.com/Baharis/hikari/commit/7400d26c6cb4304caf7d1bedccb73984aaff9a0b))

- Fix python version 3.10 being interpreted as float
  ([`e9a55d1`](https://github.com/Baharis/hikari/commit/e9a55d13a64a843541a7d2c16bba165fd988df60))

- Fix temporary completeness statistics test
  ([`e44fe89`](https://github.com/Baharis/hikari/commit/e44fe89d83d850ec5f293466ffb5fc158434c4fe))

- Fix tests of similarity index comparisons to ignore heared lines
  ([`9f0907d`](https://github.com/Baharis/hikari/commit/9f0907d4fe6fcfcf0647d48a9a123d49296fb6ac))

- Fix typo in `test_space_groups_picked_match_hall`
  ([`22594c0`](https://github.com/Baharis/hikari/commit/22594c0594fa5054c85938a8cd5718323d7d93f3))

- Fix wrong definition of d in `HallSymbol`s
  ([`0d89742`](https://github.com/Baharis/hikari/commit/0d897424401cdadd0a86f1d8cc52d914ead671f4))

- In `catalog`, change pickling protocol to 4 to support python 3.7 and 3.8 loader
  ([`8016a9c`](https://github.com/Baharis/hikari/commit/8016a9ca361e37d294c34c1428f3079faf0dd9b7))

- In `catalog`, replace `f'{var=}'` syntax from python 3.9 with `f'var={var}'`
  ([`a20a175`](https://github.com/Baharis/hikari/commit/a20a17556b06beb916a18167277adcad6a308462))

- In `resources`, replace `importlib.resources.files` from python 3.9 with `open_binary` &
  `open_text`
  ([`2c46c6c`](https://github.com/Baharis/hikari/commit/2c46c6c833199bda46280de0caf0cd75b576df2d))

- In type hints, replace `collections.abc.Buffer` from python 3.12 with simple `bytes`
  ([`0b27187`](https://github.com/Baharis/hikari/commit/0b271876b11b978c3f2b2b29bece32de1a126438))

- Increased `HklKeyRadius` precision to float64
  ([`f41abe8`](https://github.com/Baharis/hikari/commit/f41abe8c61af920e99bc7d1b889dfd689629c82b))

- Make all docstrings and regexes with backslash raw
  ([`af20f77`](https://github.com/Baharis/hikari/commit/af20f772014c8e4277039be1336b639f919b7250))

- Make all docstrings and regexes with backslash raw 2
  ([`278c5ee`](https://github.com/Baharis/hikari/commit/278c5eeadd5970c9ccfd0f304458b681e9922e11))

- Merge branch 'development'
  ([`aa8476b`](https://github.com/Baharis/hikari/commit/aa8476b194732823033da068ecabae8581233a00))

- Merge pull request #38 from Baharis/futureproof
  ([`9a96455`](https://github.com/Baharis/hikari/commit/9a964552c4bc56f945cc8565be158e30266d3a0e))

Futureproof

- Merge pull request #39 from Baharis/futureproof
  ([`87bfc47`](https://github.com/Baharis/hikari/commit/87bfc47aeca1f2d99204ca4016f30fc09e85aec1))

Futureproof

- Merge pull request #40 from Baharis/hall_symbols
  ([`cd58abe`](https://github.com/Baharis/hikari/commit/cd58abea4fa7c38fc9524a47d9bfad6ceebf5020))

Increase overall quality of handling, storing, selecting Groups; Increase user experience quality
  and code safety; Bug fixes.

- Merge pull request #41 from Baharis/development
  ([`36888c6`](https://github.com/Baharis/hikari/commit/36888c68d7ac875c55c597bb26c2c6e5f7fcfed3))

:twisted_rightwards_arrows: Bump minimum python to 3.9 and manage versioning with poetry

- Merge remote-tracking branch 'origin/futureproof' into futureproof
  ([`08778c9`](https://github.com/Baharis/hikari/commit/08778c91fc71ef2836628113f5d1925a4a572a0f))

- Minor style and docstring improvements
  ([`f7258a8`](https://github.com/Baharis/hikari/commit/f7258a87cbfe0eaf024616676aa967e99b456ffd))

- Modify implementation to account for point groups without '-p' at the beginning
  ([`12900b7`](https://github.com/Baharis/hikari/commit/12900b7cab88a2e5c79ebf86e62b04ffd94c909c))

- Most of the HallSymbols work, some errors still to fix
  ([`f6e52bc`](https://github.com/Baharis/hikari/commit/f6e52bcbd3ea9386adb097ca8a707a5314c7706f))

- One test is unpredictable on ubuntu-latest, no idea why
  ([`707b140`](https://github.com/Baharis/hikari/commit/707b1400cae5d1cea06fd7a3882e03b0fb06a344))

- Only check "older" python3.6 on ubuntu-20.04
  ([`1505106`](https://github.com/Baharis/hikari/commit/150510641dd8714c9db5142cd78fe4354251f93f))

- Preserve hkl key 'equiv' to consistently be of dtype np.int64
  ([`3bf53f6`](https://github.com/Baharis/hikari/commit/3bf53f659cf7f1cf679bff846134c1b6e3e5f4f4))

- Preserve hkl key 'equiv' to consistently be of dtype np.int64
  ([`3724333`](https://github.com/Baharis/hikari/commit/3724333f28e4ba3c7d8348994fdc7d396cc45fe3))

- Prevent deprecation failure by adding `observed=False` to hkl table groupby
  ([`fa05b53`](https://github.com/Baharis/hikari/commit/fa05b53e8fcc6fc603914873d7d1b664d6c9f2ce))

- Regenerate pickles and update pickle-dumping/loading mechanism
  ([`841e6f3`](https://github.com/Baharis/hikari/commit/841e6f3dc12c45a1ba28a1886aa8438ffa87f39d))

- Remove `HklFrame.keys` and `HklIo.keys` altogether
  ([`de8a67d`](https://github.com/Baharis/hikari/commit/de8a67dea33a65a6155c397da54a1a8ce3f25dc1))

- Remove all usage of pickle from the code
  ([`4362cdd`](https://github.com/Baharis/hikari/commit/4362cdd74220dde1f7c5c3c7d52c97db1094b33b))

- Remove debug print statements
  ([`d6d5f0f`](https://github.com/Baharis/hikari/commit/d6d5f0f065b03e180fff5ca9623bc266cc4317d8))

- Remove ReST references to deleted modules, minor fixes
  ([`6aebac0`](https://github.com/Baharis/hikari/commit/6aebac02b322207c9e6e2fb11dffaf53fcdfd1a5))

- Replace deprecated pandas `delim_whitespace` with `sep='\s+'`
  ([`5f91c8f`](https://github.com/Baharis/hikari/commit/5f91c8fea82cf6807384b329495f43430ff1ddf5))

- Rewrite `resources` to use `importlib.resources`. Drop Python 3.7 support for resources and
  `__future__.annotations`.
  ([`cf8a606`](https://github.com/Baharis/hikari/commit/cf8a606941006534067f69f22e88021074b462f7))

- Since we drop `__future__.annotations`, replace Tuple, List, OrderedDict with tuple, list, dict
  ([`c9f53b9`](https://github.com/Baharis/hikari/commit/c9f53b9733666f566700c932d891422988fcb172))

- Skip tests that rely on tk GUI for development purposes until Python 3.13.1 fixes init.tcl
  ([`b5da0ea`](https://github.com/Baharis/hikari/commit/b5da0eab41fb2d1d7efe6295cd9e3e8722b41105))

- Switch serializer from pickle to JSON to avoid version conflicts and increase security
  ([`81d927c`](https://github.com/Baharis/hikari/commit/81d927c9a36aba1b26479ba8f04a7cbc842c87ed))

- Tests fail @ Python 3.10, TODO: fix matplotlib deprecation
  ([`cabf95b`](https://github.com/Baharis/hikari/commit/cabf95bf3797656238794a47e43990ce5d7faa26))

- Todos; fix `BaseFrame` reading angles from incorrect fields
  ([`88e6380`](https://github.com/Baharis/hikari/commit/88e6380fe6cfc5c6651bea3f375661cb844b7bbb))

- Tweak, move `HallSymbol` to separate file
  ([`5ab853c`](https://github.com/Baharis/hikari/commit/5ab853ccac641d0ba68dcbc52bc7bdb4356326d5))

- Update SG pickle, simplify its dumping, add strong Group tests
  ([`2168209`](https://github.com/Baharis/hikari/commit/216820987e929967c60e121b4004a20dc5ce52d9))

- Wip: Start implementing generating SGs from Hall Symbols
  ([`2aba1de`](https://github.com/Baharis/hikari/commit/2aba1de3f9d4b90e2baa951db3b305022d0177ed))

- Write docstring, ultimately remove `point/space_group.py/.json`, adapt tests
  ([`f6fa4f0`](https://github.com/Baharis/hikari/commit/f6fa4f0df3afdf193f547b222e859330384870e7))


## v0.2.3 (2022-08-12)

### Other

- `animate_similarity_index` now accepts list of transformations instead of single step
  transformation
  ([`3574282`](https://github.com/Baharis/hikari/commit/35742821f744c9dc8d7bbd3f5b89b1f9a23e64a3))

- `animate_similarity_index` now fixes axes to same length and shows SI plot < 10 only
  ([`687cd61`](https://github.com/Baharis/hikari/commit/687cd612953ca532b65f078b7f252798879fce04))

- `animate_similarity_index` now plots using four individual axes
  ([`0303209`](https://github.com/Baharis/hikari/commit/03032090f61a46cbc20c9202ee7400a418f0ebcf))

- `hklframe.fill` now uses 4x less memory when creating initial hkl list to cut
  ([`2be4e4d`](https://github.com/Baharis/hikari/commit/2be4e4dde741ab0122179089d7bb317c4f54862a))

- Added `animate_similarity_index` for making simple animations
  ([`559d62d`](https://github.com/Baharis/hikari/commit/559d62dc8d1a5a558661fd8fb8e1e9f7ca59260a))

- Added `animate_similarity_index` to `scripts.__init__`
  ([`223f1f2`](https://github.com/Baharis/hikari/commit/223f1f2ec94521493f884b974b3549b5ab6a486f))

- Added a dedicated test case for `HklFrame.fill`
  ([`c7db318`](https://github.com/Baharis/hikari/commit/c7db3186fd7f1628fbe8b0c42c087ce1fcb7aa89))

- Added description of `m80` format to `HklIo.format` docstring
  ([`fbcbeb0`](https://github.com/Baharis/hikari/commit/fbcbeb0a1a32880e74bf841cd5221e9d425c41aa))

- Added docstring and fixed text placement mistake in `animate_similarity_index`
  ([`704e65b`](https://github.com/Baharis/hikari/commit/704e65b0782735dd2e8aae0f042512be3e789a80))

- Bumped hikari version to 0.2.3
  ([`304ee6b`](https://github.com/Baharis/hikari/commit/304ee6b599fc5f60b53e3b48934b59f9328a5b71))

- Changed try/except loops to catch errors when reading hkl whenever fixed-width line is too short
  ([`c28f1d5`](https://github.com/Baharis/hikari/commit/c28f1d50ba6eb5a9a9b544e0cb8fc870328fe705))

- Created `utility.numpy_tools.str2array` to read compact matrix representations
  ([`3547788`](https://github.com/Baharis/hikari/commit/3547788bbf7f4f01b1c326cb54e66153fbe910b6))

- Created test suite for `utility.numpy_tools.str2array`
  ([`5762401`](https://github.com/Baharis/hikari/commit/5762401e81d1d804ae4a74faf477cb15818ce9a6))

- Defined `HklFrame.fill2` method in order to optimise memory use when filling with reflections
  ([`b957ff8`](https://github.com/Baharis/hikari/commit/b957ff8800f98118961ee17ea84d2d6d8971b722))

- Defined a basic structure of the `Selling` class
  ([`2e13f9c`](https://github.com/Baharis/hikari/commit/2e13f9c38b4fc3d8feb38d7abb11f8d12747873f))

- Defined a very bare-bones m80 reflection file type in `hkl_formats_defined.json`
  ([`da566de`](https://github.com/Baharis/hikari/commit/da566de1d9baf0985d527b6cf9de76cc81992f91))

- Defined an empty dummy key type `_` in `HklKeys` to hold irrelevant data
  ([`4fd12dc`](https://github.com/Baharis/hikari/commit/4fd12dcc7547281792758cae8e28bb4d4f78282c))

- Fixed and tested `HklFrame.fill2`. Now it won't produce hkl > HKL_LIMIT and is 2-3x faster
  ([`628194b`](https://github.com/Baharis/hikari/commit/628194b7daf44c6795002da69b289a632174112d))

- Implemented `HklFrame.fill3`, which is still 2-3x faster and severely less memory hungry than
  `fill2`
  ([`be9a05a`](https://github.com/Baharis/hikari/commit/be9a05a5947b54b92e9ea52ee5a224c8e45189a1))

- In `dataframes.base` defined new class `Selling` and its matrices to handle Selling reduction
  ([`cc3ca38`](https://github.com/Baharis/hikari/commit/cc3ca38c333bc0215b578901fc5ebdbdfd3f47d0))

- In `HklFrame.fill`, in `hkl_walls` changed `l` to `l_` to avoid ambiguity
  ([`89d25ca`](https://github.com/Baharis/hikari/commit/89d25ca62d5fca41707b4cdf3a8cb81aa1b4225f))

- Increased `HKL_LIMIT` to 127. New fill uses less memory and is 9x faster for very oblique cells
  ([`d5f02f4`](https://github.com/Baharis/hikari/commit/d5f02f4467fd665dc5b97b28bb6ccd931c51bded))

- Increased readability of `HklFrame.fill`
  ([`e5e0411`](https://github.com/Baharis/hikari/commit/e5e04115196f266e6e5905fd679e5fa83d03782a))

- Merge pull request #34 from Baharis/animation
  ([`f19ef43`](https://github.com/Baharis/hikari/commit/f19ef435ec0ef40e185d982c6e1b45165cdc5c2c))

Animation

- Merge pull request #35 from Baharis/m80
  ([`0053f2d`](https://github.com/Baharis/hikari/commit/0053f2d5a4eb52286fe1bea651b8ad556f3c5326))

M80

- Merge pull request #36 from Baharis/animation
  ([`0e12e30`](https://github.com/Baharis/hikari/commit/0e12e309e14f1fea94432df7f74e573f26f32297))

Similarity index animations fix

- Merge pull request #37 from Baharis/reduction
  ([`5f0fa00`](https://github.com/Baharis/hikari/commit/5f0fa0098d8af670898a04a5bc184380cd4756f3))

Reduction

- Moved Selling reduction tools from `dataframes.base` file into `BaseFrame`
  ([`2ab2b2c`](https://github.com/Baharis/hikari/commit/2ab2b2c46692b294ccfd4cc9d2062456d16dd187))

- Redefined empty/trash entry type in hkl files from `_` to `None`
  ([`f8eda14`](https://github.com/Baharis/hikari/commit/f8eda141dc0cae9a7b823b6d60ddfbb7d1062092))

- Removed attempts to define partially reduced cell - maybe useful in the future, too much work now
  ([`c6282ac`](https://github.com/Baharis/hikari/commit/c6282ac9af0ac13643ea5a6ed70f1efa9fdf2d25))

- Renamed `reduction` to `semi-reduction`, as it doesn;t need unique and fixed ordering
  ([`9cc6de4`](https://github.com/Baharis/hikari/commit/9cc6de406f8fab6ecd29e1bb4f402ee0afe6a255))

- Started animate_similarity_index based on so #7819498 question
  ([`526c7fc`](https://github.com/Baharis/hikari/commit/526c7fca56730c81c3d8ec025e43f9207e6d8a89))

- Substituted old `fill` with new, faster and less demanding one
  ([`9fb1b6f`](https://github.com/Baharis/hikari/commit/9fb1b6fa8fe8b636f57c1e1744480f01817c328b))

- Tested `HklFrame.fill3` tu be 2-3x/~1.5x faster in problematic/standard case, Replaced `fill`
  ([`3e042af`](https://github.com/Baharis/hikari/commit/3e042af6d1dd25014e8fdf8b1f6e8aebfc08e0ff))


## v0.2.2 (2022-06-20)

### Other

- `angular_explorer` now prints values of property in focus points to the listing file
  ([`501780f`](https://github.com/Baharis/hikari/commit/501780f05e3006b48270fc3877e31ed998e34c53))

- `baseframe.fill_from_cif_block` now correctly raises KeyError on missing data if fragile
  ([`16f0fcf`](https://github.com/Baharis/hikari/commit/16f0fcfa0bb45c5a423c997760fc5d7549b38d17))

- `cifblock.get_as_type` now correctly returns unmodified default if given
  ([`09ae41d`](https://github.com/Baharis/hikari/commit/09ae41d3c733dd05c6cdb9808804e30b93dedb12))

- `cifframe.get_as_type` now correctly does not raise KeyErrors
  ([`63445f8`](https://github.com/Baharis/hikari/commit/63445f8200ab621fdaae611096a696c371583c88))

- `cifreader.databuffer` moved to `CifIO` and made abstract
  ([`6c8af86`](https://github.com/Baharis/hikari/commit/6c8af86be9b1015f70a5e0edcc44b79543bcb62d))

- `cifreader` now does not store last '\n' in multilines; accordingly, `CifWriter` adds it
  ([`110ccaf`](https://github.com/Baharis/hikari/commit/110ccaf958407059465f46d43a4444a0527a867f))

- `cifwritebuffer` now correctly clears lists when flushing
  ([`b854b51`](https://github.com/Baharis/hikari/commit/b854b510916b8b3da88d5bbe8792281d5037955a))

- `cifwriter.write` method made to work after some tests
  ([`9b3c497`](https://github.com/Baharis/hikari/commit/9b3c497dcfc0a70ef1fc53ba191b84b6ab4ba4e8))

- `cifwriter.write` method now writes the last data of each block correctly
  ([`fd43fd9`](https://github.com/Baharis/hikari/commit/fd43fd9de20538808f545366582d8a7dce48c8a1))

- `cifwriter.write` now forcibly enquotes single multiline cif entries
  ([`b49d292`](https://github.com/Baharis/hikari/commit/b49d292742bfaad0fd39f297356ab2bc9b6876fb))

- `cifwriterbuffer` now clumps cif-core-missing entries together if they are similar enough
  ([`3c6a90a`](https://github.com/Baharis/hikari/commit/3c6a90aa81d07abda97ec3270c5959fee55bf0f8))

- `hikari.scripts.calculate_similarity_indices` now returns 1 for negative adps
  ([`41c8e93`](https://github.com/Baharis/hikari/commit/41c8e93bac76e740a57b828d9adced4dbeb13ed8))

- Added `MULTILINE_QUOTE_REGEX` to prevent commented "\ndata" from triggering as block
  ([`acaf5c9`](https://github.com/Baharis/hikari/commit/acaf5c9651a2c6cf5b8bb85a66db43337da90064))

- Added a `NaCl_negative_ADPs` block to test NaCl cif file
  ([`8694104`](https://github.com/Baharis/hikari/commit/8694104ca1acce3a1709622df9af802bb675401d))

- Added a bunch of tests for `hikari.dataframes.cif` objects
  ([`808f72b`](https://github.com/Baharis/hikari/commit/808f72b229343c1197903ceea1197420fe9b6e62))

- Added a few tests for `hikari.scripts.calculate_similarity_indices`
  ([`27e21b2`](https://github.com/Baharis/hikari/commit/27e21b204a0fe8c97b0046f00a76fa5aeb5ec274))

- Added to-be-tested `write` methods to `CifFrame` and `CifBlock`
  ([`2f19333`](https://github.com/Baharis/hikari/commit/2f19333dfc0c279fe2a7076e219c697e1df6d7b6))

- Added two simple tests for creating `potency_map`s
  ([`432147d`](https://github.com/Baharis/hikari/commit/432147d2311155e4429cfd758153f5257086432a))

- Added two simple tests for writing cif files in a reproducible way
  ([`f7de3cb`](https://github.com/Baharis/hikari/commit/f7de3cb6d65e46ef9517431b9dd390ee5f11d8d7))

- Bumped `hikari` version to 0.2.2
  ([`5796393`](https://github.com/Baharis/hikari/commit/5796393c10116c89187656b6e7606a7cffab14aa))

- Cleared TODOs and misleading if-name-main
  ([`04088e8`](https://github.com/Baharis/hikari/commit/04088e8c0b4f4e991ff41c9bc4dc911bd9f5e81b))

- Fixed `CifReader` and `CifWriter` docstring to generate documentation correctly
  ([`c659e45`](https://github.com/Baharis/hikari/commit/c659e452c94d7265b85f61cf67f28d26ce1d0da7))

- Fixed error with quoted comments by protecting quotes before stripping comments
  ([`3c139ba`](https://github.com/Baharis/hikari/commit/3c139ba83fb47e495d3dc678582e2e82c7dcf964))

- Fixed looped multilines being incorrectly assigned to loop_keys in `CifReader`
  ([`688ef19`](https://github.com/Baharis/hikari/commit/688ef199f22dcd7c5126992e63c5df72b9b20be4))

- Fixed multilines being one newline too long after `CifWriter`
  ([`39048bd`](https://github.com/Baharis/hikari/commit/39048bd92e497f03a51e9c301a46134d3eea950c))

- Implemented `CifValidator.get__category` to check groups in cif validator
  ([`4434883`](https://github.com/Baharis/hikari/commit/4434883ed2c36b690b6ace9f0b203ac0a89b30a4))

- Implemented `CifValidator.get__list` to check if `_list` is 'yes' in cif dict
  ([`09c199f`](https://github.com/Baharis/hikari/commit/09c199f22992c2ef226da03a44cad76837c729b8))

- Implemented `CifWriterBuffer` to write individual data elements to cif file
  ([`ebd2a3f`](https://github.com/Baharis/hikari/commit/ebd2a3fc550ac5d4ee1d6b6028270b662538a111))

- Merge pull request #31 from Baharis/CifIO
  ([`0fe1d1a`](https://github.com/Baharis/hikari/commit/0fe1d1a41f1badc065cad2e6293b2f4d6fc14cfc))

Cif io

- Merge pull request #32 from Baharis/similarity
  ([`97056be`](https://github.com/Baharis/hikari/commit/97056be67686bad32c7446f242a2ef6bf91b427c))

Similarity

- Merge pull request #33 from Baharis/property_at_focus
  ([`2cbd25c`](https://github.com/Baharis/hikari/commit/2cbd25c42bcc5fb8d96e8c76cfda1d5234f3bce4))

Property at focus

- Modified regex expressions used in `CifReader`
  ([`a9e0ebc`](https://github.com/Baharis/hikari/commit/a9e0ebc8d6267f44a5e93eac99d3a64f52185732))

- One can now get hidden keys from `CifValidator` (such as `_atom_site_aniso_B_11`) without
  expanding
  ([`4a2368c`](https://github.com/Baharis/hikari/commit/4a2368cb6cf41e60ea08699d2609b02412fc6385))

- Overloaded `contains` method of `CifValidator` based on results of `get`
  ([`cfe5935`](https://github.com/Baharis/hikari/commit/cfe593562ae3e8b9e83d39969cb40e7fa3e7815b))

- Registered "\n" as whitespace and added `DATA_BLOCK_REGEX` to CifReader
  ([`80a2663`](https://github.com/Baharis/hikari/commit/80a26638787d5d35f7b5ca816ffdd55b1a2bd5e6))

- Removed "# noinspection PyTypeChecker" from code to prevent type misinterpretation by git
  ([`8f06ad4`](https://github.com/Baharis/hikari/commit/8f06ad40d52e63d946df632f9e3b17360589da5d))

- Restricted `CifReader.DataBuffer` (now `IOBuffer`) to a very general abstract class
  ([`39beec0`](https://github.com/Baharis/hikari/commit/39beec0b2e724740c2e20c096f066de3b4c4a35e))

- Restructured buffer classes out of IO classes; they are going to play a large role
  ([`9e68088`](https://github.com/Baharis/hikari/commit/9e68088573c41b86091511c45e8a94426f800d1f))

- Shortened regex expressions, found error with quoted pound signs
  ([`24447c8`](https://github.com/Baharis/hikari/commit/24447c8a51c3c8776e1ebf596fc80d7919955b76))

- Simplified `calculate_similarity_indices` docstring using type hints
  ([`67a0834`](https://github.com/Baharis/hikari/commit/67a0834f9586c550bfaad107de2683c4900c8a04))

- Simplified `CifReader` and `CifReaderBuffer` by protecting multilines
  ([`770c79a`](https://github.com/Baharis/hikari/commit/770c79afd677abe23ac0041d4032a7f737ebe15a))

- Since `CifBlock.get_as_type` is now allowed to return `None`, modify tests to check that
  ([`ff172ec`](https://github.com/Baharis/hikari/commit/ff172ecdf4ad0ac17f65c5ae47f8dc11db5037f1))

- Slightly optimised the regex
  ([`0f8eab7`](https://github.com/Baharis/hikari/commit/0f8eab79c43a03a9deaf49c2aa57de798aabb12b))

- Started creating `CifWriter` and `CifWriter.WriterBuffer`
  ([`f0fdb52`](https://github.com/Baharis/hikari/commit/f0fdb52359fe2d4329251a0f99d1251f4a2d9e39))

- Started reconstructing `CifReader` using more regex expressions
  ([`dc7a471`](https://github.com/Baharis/hikari/commit/dc7a4716b18a9f5f788ac1c53bcbde9d65466cdf))

- Successfully tested `calculate_similarity_indices` script against MW's datasets
  ([`4b52e5f`](https://github.com/Baharis/hikari/commit/4b52e5fd3e2d368de4153a3c13b39e07641a7e3f))


## v0.2.1 (2022-06-10)

### Other

- Added missing package data to the setup file
  ([`b81726d`](https://github.com/Baharis/hikari/commit/b81726d3d5cfcbe7ace75b420166f8a7513f3b5b))

- Bumped `hikari` version to 0.2.1
  ([`a1f868b`](https://github.com/Baharis/hikari/commit/a1f868be97511e7b2cb75fae2d1a0f85a20f9a20))

- Merge pull request #30 from Baharis/package_optimisation
  ([`e460d25`](https://github.com/Baharis/hikari/commit/e460d2577e63920410343846e01f5b248eedca29))

Resources hotfix


## v0.2.0 (2022-06-09)

### Other

- 'cifblock.get_as_type' now correctly turns lists into lists instead of maps
  ([`45ce7b2`](https://github.com/Baharis/hikari/commit/45ce7b2901101cd1ac71f64cb9b33629777fd8d1))

- 'hikari.dataframes.baseframe' can again read cif blocks using `fill_from_cif_block`
  ([`8c025b2`](https://github.com/Baharis/hikari/commit/8c025b2d2c0831f94c3a37040a36d0cfbeb2fd62))

- `baseframe.imported_from_cif` is now name-indexed dict of 3-els instead of list of 4-els
  ([`f3a8b09`](https://github.com/Baharis/hikari/commit/f3a8b096f44452304eb8aa908eb50812500e14b7))

- `cifframe.get_as_type` uses default value now, raises KeyError if no default nor value found
  ([`be25b58`](https://github.com/Baharis/hikari/commit/be25b5829b02198455f079e484aba74e47bba30f))

- `cifio.__init__()` now checks whether item should be a list based on length and `CifValidator`
  data
  ([`702bd5d`](https://github.com/Baharis/hikari/commit/702bd5d98763456dec844ba9b0141159121d731e))

- `cifio` now correctly reads singletons as lists via `CifValidator`
  ([`6c02c7e`](https://github.com/Baharis/hikari/commit/6c02c7e488c1bd696bddd440f5464d71eb95a387))

- `cifreader` now correctly ignores comments by calling `strip_comments` when reading
  ([`788823f`](https://github.com/Baharis/hikari/commit/788823f8958072135852e63b665fbec9ef267123))

- `det3x3` renamed from `udet`, moved from `UBaseFrame` to `hikari.utility.mathtools`, and made used
  by `BaseFrame`
  ([`0fd3109`](https://github.com/Baharis/hikari/commit/0fd3109db1cc77441c5d756d2df16d8235e77d03))

- `hikari.dataframes.cif.cifblock` now properly raises Error when no default nor value are found
  ([`cb9fa7a`](https://github.com/Baharis/hikari/commit/cb9fa7a4d564516f62ac628ccc939672697eadc7))

- `revert_whitespace` now acts on single string; Created `CifWrite` from CifIO
  ([`68348a7`](https://github.com/Baharis/hikari/commit/68348a7644e5f24d8b219a45e16d6f58ce0d2ed1))

- `test_dataframes` constants now are correctly uppercase
  ([`9d1405a`](https://github.com/Baharis/hikari/commit/9d1405ab9ded55e39ee811b49a82c9dfdc530663))

- `ubaseframe` tested, moved to its own file, and added to `hikari.dataframes.__init__`
  ([`5ef6a7d`](https://github.com/Baharis/hikari/commit/5ef6a7d2987a1be212990eca77eaea09bdc9896d))

- `utilities.deg2rad` can now /correctly/ handle ufloats
  ([`ad9b880`](https://github.com/Baharis/hikari/commit/ad9b880344523c5e6311e926f7c16ffe01f5f35e))

- `utilities.deg2rad` can now handle ufloats
  ([`8e91472`](https://github.com/Baharis/hikari/commit/8e9147277595db1f88ca5c984edecb456dd8bad8))

- `utility.det3x3` now uses indexing instead of slicing to properly communicate with `unumpy.matrix`
  ([`478d111`](https://github.com/Baharis/hikari/commit/478d11182c32d08ad224c182f02127680d34b64e))

- Adapted `README.md` to include e.g. fcf/cif support
  ([`412e8ae`](https://github.com/Baharis/hikari/commit/412e8aec2c50f7840e0e98addf5326c690920960))

- Added `__init__.py` to test folder for coverage integration
  ([`c5a350f`](https://github.com/Baharis/hikari/commit/c5a350f760b54f75f0f458d696e9a7c1e05ba327))

- Added `CifBlock` to `hikari.dataframes.__init__` to accessible from outside
  ([`26b739b`](https://github.com/Baharis/hikari/commit/26b739b5c31767e12ac5cc65b87a3171c732674e))

- Added `compare_adps` to `hikari.scripts` submodule
  ([`139826e`](https://github.com/Baharis/hikari/commit/139826e7b2cea3964c56d54efafce4f2cf2edd6e))

- Added `get_as_type` method to `CifBlock`
  ([`b5f8a3f`](https://github.com/Baharis/hikari/commit/b5f8a3f7a73d65b4bcf6ed6f28c46541c6c4e9eb))

- Added `nacl_cif`, `nacl_fcf`, and `cif_core_dict` to `hikari.resources`
  ([`edfbede`](https://github.com/Baharis/hikari/commit/edfbede5d102a6c947399c593e8f0f2d2771925e))

- Added `read` method to `CifBlock`
  ([`c74a1eb`](https://github.com/Baharis/hikari/commit/c74a1eb96d1de53d847c4c62b82cf04c2c74748c))

- Added a `CifValidator` object to `hikari.dataframes.cif`
  ([`c3a5e4a`](https://github.com/Baharis/hikari/commit/c3a5e4a063834de29972dca51f3f4fbca07bf7ca))

- Added a `common_prefix` function to `hikari.dataframes.cif`. Cif reader needs rewriting...
  ([`4fe2cce`](https://github.com/Baharis/hikari/commit/4fe2ccea475538ed2a89df26ca5ac67d6fea42b8))

- Added a generic script for calculated similarity index for one atom
  ([`be8341a`](https://github.com/Baharis/hikari/commit/be8341a47a0bb19d6ac1dedb97f9cb09f0a55e9b))

- Added a shameful codecov badge
  ([`7735cca`](https://github.com/Baharis/hikari/commit/7735ccac893ec82e1fdca898996434b451000c4a))

- Added a simple `reformat_hkl` script to match readme
  ([`aed1335`](https://github.com/Baharis/hikari/commit/aed1335863f3e4b179770c634fb2ad8385b67942))

- Added and commented some basic script for many-atom `compare_adp` script
  ([`bcb9b57`](https://github.com/Baharis/hikari/commit/bcb9b57dde09f7591abdfe43b2bddaa64b101ddf))

- Added and updated tests for `BaseFrame`, `CifBlock`, `CifFrame`, and `UBaseFrame`
  ([`03aa1ad`](https://github.com/Baharis/hikari/commit/03aa1ad14886bd082ac5b3addb6a5b5a750daa59))

- Added CI tests badge
  ([`ce8f2c7`](https://github.com/Baharis/hikari/commit/ce8f2c7d3c2ff9ce12b41db635945633e2a154e5))

- Added cif dictionary to the resources
  ([`4b0f21c`](https://github.com/Baharis/hikari/commit/4b0f21c2b72495eb6b557b3f9cd5445119a9bb51))

- Added codecov workflow
  ([`c5f1ec3`](https://github.com/Baharis/hikari/commit/c5f1ec320642211154cf534b9b5963a3271ef55b))

- Added docstrings to `hikari.dataframes.CifIO` methods
  ([`83756a3`](https://github.com/Baharis/hikari/commit/83756a3c1faf76a53b2bdf16bc5d9a82f9313161))

- Added inactive `validate` to `CifIO.__init__()` to prevent `CifValidator` infinite loop
  ([`482d69b`](https://github.com/Baharis/hikari/commit/482d69bf861863a692deab307760cc247e9d2964))

- Added other NaCl files to the resources for the purpose of testing
  ([`67b0d16`](https://github.com/Baharis/hikari/commit/67b0d1683e121ed4b7a84757bdef452ca5d41689))

- Added project `requirements.txt` to Sphinx requirements
  ([`270df0a`](https://github.com/Baharis/hikari/commit/270df0ac89b9037a2cd447c35f2449c0b995d8a7))

- Added readthedocs files to try and fix empty documentation problem
  ([`534ab09`](https://github.com/Baharis/hikari/commit/534ab0914df0626834c7e31727e76d07bf234b02))

- Added reminder tofix the documentation later on
  ([`b09e5c0`](https://github.com/Baharis/hikari/commit/b09e5c0b0cff098abf06458412cabafc9f9b658f))

- Added requirement for `uncertainties` package version 3+
  ([`550e1fe`](https://github.com/Baharis/hikari/commit/550e1fe8181f99d3d115947fbe81c130960a4219))

- Added tests for hikari.utility.certain_float submodule
  ([`7cdd941`](https://github.com/Baharis/hikari/commit/7cdd94162d7af97eeb9b4d92b9e0474d54bc39f5))

- Added, updated, corrected tests for `CifBlock`
  ([`5efb750`](https://github.com/Baharis/hikari/commit/5efb7505bda9e29f4f030b0b8a4b82fd3dd0b298))

- All cell parameters of `UBaseFrame` are now ufloats
  ([`6511006`](https://github.com/Baharis/hikari/commit/6511006ef8944896e8adfd0c5d9b5526be4b7980))

- Bumped `hikari` version to 0.2.0
  ([`ff31b56`](https://github.com/Baharis/hikari/commit/ff31b568974207c20e35f808ea35c018b1988c58))

- Bumped the version to check whether this impacts documentation
  ([`301a492`](https://github.com/Baharis/hikari/commit/301a4923b62f76b812ef7c0a86fb684135d3996d))

- Changed `CifFrame` class to be a subclass of `OrderedDict`
  ([`d5c42b8`](https://github.com/Baharis/hikari/commit/d5c42b819b0e85cc6f224491c540fd8ff9801c50))

- Changed `CifReader.MATCHING_QUOTES_REGEX` to disregard quotes not followed by whitespace
  ([`a4b8dbc`](https://github.com/Baharis/hikari/commit/a4b8dbc13e40f17fded1515aed0ab36685e2b886))

- Changed `CifReader.MATCHING_QUOTES_REGEX` to not count `''` as quote
  ([`9ec480d`](https://github.com/Baharis/hikari/commit/9ec480d1e0463a5b22576c053535d6830e0b74c2))

- Changed `CifReader` to be able to read cif dictionary, but with problems
  ([`78f65ce`](https://github.com/Baharis/hikari/commit/78f65ce9731e698eb43e74b85271d225796feb1d))

- Changed `CifReader` to correctly `protect` and `release` quoted words
  ([`e654168`](https://github.com/Baharis/hikari/commit/e654168e278004f69f0ce7115e962ff2ee427c3b))

- Changed `test_dataframes.py` so that files are created only once
  ([`55dabc7`](https://github.com/Baharis/hikari/commit/55dabc74024edbf5afd3ff0f39bd41b7ca7d61a1))

- Cleaned `BaseFrame` and `UBaseFrame` inits' code
  ([`7244d2f`](https://github.com/Baharis/hikari/commit/7244d2f2d25f14259dc75d765adfdf499e7cf141))

- Corrected grammar in `README.md`
  ([`c2a1ee3`](https://github.com/Baharis/hikari/commit/c2a1ee3e25bf0fbb8352d9c707ea0453f1ad429d))

- Corrected order inside `CifReader.substitute_quoted_whitespace`
  ([`22ce28d`](https://github.com/Baharis/hikari/commit/22ce28d20104def681f5ba594a627377c46483cf))

- Create .readthedocs.yaml
  ([`ee0a18c`](https://github.com/Baharis/hikari/commit/ee0a18c641b9a0ce484a63b1169fdf01744efac3))

- Created empty `CifBlock` class to handle individual blocks of `hikari.dataframes.CifFrame`
  ([`c00ef01`](https://github.com/Baharis/hikari/commit/c00ef019d365b656406ea5fb1d429dc0e41767d8))

- Defined `cfloat` "type" (str2ufloat2float converter) to `hikari.utility`
  ([`edd2f51`](https://github.com/Baharis/hikari/commit/edd2f518c3bcafa7f68c21d7dfcca6fd76b2ed7a))

- Documentation corrected throughout the project to include empty line between descriptions and
  parameters
  ([`48b57dc`](https://github.com/Baharis/hikari/commit/48b57dc18d9694e4883a099dc066f198ef3f4545))

- Extended tests for `hikari.dataframes` submodule
  ([`5bd29d2`](https://github.com/Baharis/hikari/commit/5bd29d27e71b501142e3790796226d4ee89ee2b1))

- Extended tests for `hikari.utility.math_tools` submodule, changed bug in `rotation_from(a, to=-a)`
  ([`a4c29d3`](https://github.com/Baharis/hikari/commit/a4c29d3df9069abcea45ef8e12e99c9f337f6deb))

- Extended tests for hikari.utility.chem_tools submodule
  ([`21d77df`](https://github.com/Baharis/hikari/commit/21d77dfaa67660581f54f800fbf64b2339c23cd0))

- Extended tests for hikari.utility.interval submodule (which itself was simplified a bit)
  ([`0e3ef97`](https://github.com/Baharis/hikari/commit/0e3ef9709d59667220135e40a380e2fafb6c6cbc))

- Extracted `CifReader` out of `CifIO` class
  ([`7dfa2d3`](https://github.com/Baharis/hikari/commit/7dfa2d3bf4d7bbad70c1ef2d60bc470f2e0de006))

- Fixed `cfloat` to really return nominal value instead of ufloat
  ([`47219cf`](https://github.com/Baharis/hikari/commit/47219cfe98d13c8adc0315c245fbe6ac1584ecd5))

- Fixed `CifReader.parse` to flush `DataBuffer` at the start of a "loop_"
  ([`aa18bef`](https://github.com/Baharis/hikari/commit/aa18bef14542eee358df0c8dce6d9fdec3aa611e))

- Fixed minor formatting errors
  ([`9d0ca06`](https://github.com/Baharis/hikari/commit/9d0ca06217940fa341111fa90a60921fdf9c88e3))

- Fixed problem with premature termination of reading in `CifIO`
  ([`ab6e123`](https://github.com/Baharis/hikari/commit/ab6e123553c8b85b2cb773adac8666b9ff87cb70))

- Fixed regex and loop bugs in `substitute_whitespace_in_quotes` in `CifIO`
  ([`bb798a2`](https://github.com/Baharis/hikari/commit/bb798a23a276f45bed359627f35b16a5bf196e9f))

- Fixed unnecessary hyphen in `.readthedocs.yaml`
  ([`a5ddebd`](https://github.com/Baharis/hikari/commit/a5ddebd3489ba87cdb30e560076adb5f66a5e473))

- Implemented `blocks` property in `CifIO`
  ([`61e8e15`](https://github.com/Baharis/hikari/commit/61e8e157ad0a97bf1f2768c9880b7caf622a83c7))

- Implemented first version of `UBaseFrame`, subclass of `BaseFrame` capable of handling uncertainty
  ([`747ac01`](https://github.com/Baharis/hikari/commit/747ac01c296dc6b6dff33a87a0059ebbf531f654))

- Improved code clarity regarding floats and ufloats in `calculate_similarity_indices`
  ([`ef371f6`](https://github.com/Baharis/hikari/commit/ef371f632de3c2443c9274ab7f757623f0fff533))

- Improved output clarity in `calculate_similarity_indices`
  ([`f919218`](https://github.com/Baharis/hikari/commit/f919218cb70fe8b38817cb718c41a8ac370893bb))

- In `BaseFrame` and `UBaseFrame`, removed padding to increase maintainability
  ([`232d9a1`](https://github.com/Baharis/hikari/commit/232d9a1dea5739ca8c589a20619d53f9d47947ca))

- Included recently added submodules to documentation
  ([`04185a9`](https://github.com/Baharis/hikari/commit/04185a9ced87066c21f3680532857eeff60c3c10))

- Increased code maintainability as suggested by 'codefactor' in multiple files
  ([`80ceb53`](https://github.com/Baharis/hikari/commit/80ceb53bbc3f9982a599d1d564a1c7b3067e508a))

- Many `BaseFrame` attributes are now protected instead of private to allow changes in descendants
  ([`090e080`](https://github.com/Baharis/hikari/commit/090e080e315c5a47aafd03526727999e3637e92f))

- Merge pull request #19 from Baharis/Readthedocs-integration
  ([`69ba4ce`](https://github.com/Baharis/hikari/commit/69ba4ce4181a6c28e7464b17ec2bba32a37a5686))

Integrated repository with [readthedocs](https://hikari.readthedocs.io/en/stable/) and accordingly
  added badges to the `README.md` file

- Merge pull request #20 from Baharis/cif_frame
  ([`ef97b36`](https://github.com/Baharis/hikari/commit/ef97b36ffefdd58a2faf6027567ceab7efb2f5d3))

Cif frame

- Merge pull request #21 from Baharis/cif_frame
  ([`8a614d6`](https://github.com/Baharis/hikari/commit/8a614d6ad97be14d7941c14ba0f668f6ae5273b5))

Cif frame

- Merge pull request #22 from Baharis/UBaseFrame
  ([`294e9f0`](https://github.com/Baharis/hikari/commit/294e9f03aeadb1962808b2cd4e9aa0ce20aa3f49))

UBaseFrame

- Merge pull request #23 from Baharis/similarity
  ([`8f2e525`](https://github.com/Baharis/hikari/commit/8f2e5252db6e85048544e332b7f7a70f4d6ee56d))

Fixed the problem with requirements not being correctly read after specifying readthedocs
  configuration

- Merge pull request #24 from Baharis/similarity
  ([`4a436a9`](https://github.com/Baharis/hikari/commit/4a436a9f284b9319ab1f0aa4d902ad8561c7f833))

Similarity

- Merge pull request #25 from Baharis/similarity
  ([`2d67538`](https://github.com/Baharis/hikari/commit/2d675385a1bf814af985f079371d1f61887b91a4))

Similarity

- Merge pull request #26 from Baharis/tests
  ([`3e82126`](https://github.com/Baharis/hikari/commit/3e82126e0b426898743fd9869e7a958bfc24ee6e))

Tests

- Merge pull request #27 from Baharis/tests
  ([`561cf78`](https://github.com/Baharis/hikari/commit/561cf786d954ff3cab206fa08976638bcdac90bb))

Modified tests to accommodate codecov integration

- Merge pull request #28 from Baharis/cif_dictionary
  ([`2a184f8`](https://github.com/Baharis/hikari/commit/2a184f83672ad52604e25d5767a2eccda13641f8))

Cif dictionary

- Merge pull request #29 from Baharis/package_optimisation
  ([`00697c5`](https://github.com/Baharis/hikari/commit/00697c5281348544b3e722de27e4fc83cb0dacd0))

Package optimisation

- Merge remote-tracking branch 'origin/tests' into tests
  ([`e50d355`](https://github.com/Baharis/hikari/commit/e50d355e4fb8cd1f9a32f76a2716ed50d17a9ec0))

- Minor style fixes in `hikari.dataframes`
  ([`689c08e`](https://github.com/Baharis/hikari/commit/689c08e304d84358921229a4f11112b8155ccd74))

- Mixed too short lines in documentation
  ([`f58ff33`](https://github.com/Baharis/hikari/commit/f58ff33e5f4ddf765c51fbd3fe7cc83fcfccb380))

- Moch-import uncertainties to check whether this package does cause problems
  ([`8e4862e`](https://github.com/Baharis/hikari/commit/8e4862e106959f976f7e2a75d71ad7e59fb062fa))

- Modified `CifIO` and `CifFrame` to read all blocks instead of one specified
  ([`d904462`](https://github.com/Baharis/hikari/commit/d90446254e3b6a8335236ed2c4bdd7fa8c4351ef))

- Modified `hikari.dataframes.CifFrame` to read using new `CifIO` and added docstrings
  ([`c920e40`](https://github.com/Baharis/hikari/commit/c920e404ef7355edeec42346c8e1c2f1c5f4c0ff))

- Modified sphinx configuration files to avoid documentation errors in the future
  ([`0869493`](https://github.com/Baharis/hikari/commit/08694936029e4bf6a07abd37b51d904d6c3fdd35))

- Modified the `compare_adps` script to match its docstring
  ([`b8b23ca`](https://github.com/Baharis/hikari/commit/b8b23caffb78b8c0cd31351ccbc8d5d8e949d23b))

- Moved NaCl files to test, as they are not an essential part of the package
  ([`4f4d25d`](https://github.com/Baharis/hikari/commit/4f4d25db8cb2359de02ba37bfc9c244fbef0c09f))

- Parameter `normalize` of `calculate_similarity_index` now correctly normalises
  ([`359ec38`](https://github.com/Baharis/hikari/commit/359ec38a3d3a74052510d365d441180984a4d3be))

- Parameter `prepend` added to `CifReader.read()` to get correct block names in `CifValidator`
  ([`ae80b9f`](https://github.com/Baharis/hikari/commit/ae80b9f65bb72197e0061527e4080f6804c1c8a4))

- Prepared a docstring for `hikari.scripts.compare_adps.compare_adps`
  ([`f05d8aa`](https://github.com/Baharis/hikari/commit/f05d8aac5075cf9fbbae285213907265bbd024fb))

- Prepared a rough draft of a `CifIO` class for reading cif files
  ([`8a88b58`](https://github.com/Baharis/hikari/commit/8a88b58b891da94a0c9b19dd9a5338384371ff84))

- Prepared a successful test version of `hikari.scripts.compare_adps.compare_adps` script
  ([`bc0e655`](https://github.com/Baharis/hikari/commit/bc0e6558328059301488ebaa590920a6cc160d0a))

- Re-added doc html static path due to new docs problems
  ([`6367384`](https://github.com/Baharis/hikari/commit/6367384b6fe4aa1bc0bd573fbff7bd0de656c755))

- Re-added doc requirements and specified them
  ([`f748c1c`](https://github.com/Baharis/hikari/commit/f748c1c0364970b62c2b8438373d9806aa9ace33))

- Re-added fail on warning
  ([`2896498`](https://github.com/Baharis/hikari/commit/28964989c47a1e8a21aa35868bf3def9c92f6d98))

- Reduced `CifBlock` to a pure OrderedDict clone and made it item type in `CifFrame`
  ([`c5f767d`](https://github.com/Baharis/hikari/commit/c5f767d447ebb658cd1f72cd0c79252f77dbafab))

- Remove python requirements for building the docs
  ([`f6f2c7e`](https://github.com/Baharis/hikari/commit/f6f2c7e03368d8d5b0282dc6a726eb634d971ba5))

- Removed '_static' from sphinx' html_static_path
  ([`2860c09`](https://github.com/Baharis/hikari/commit/2860c0930564f33a7145bb106bd5232cc184f5fc))

- Removed _static reference, as it does not look like its connected to LaTeX errors
  ([`e5e0d86`](https://github.com/Baharis/hikari/commit/e5e0d86f91df1755d8f5132ba8651e821af2d1bb))

- Removed `common_prefix` function from `hikari.dataframes.cif` as its functionality was found
  dispensable
  ([`e7aa7fd`](https://github.com/Baharis/hikari/commit/e7aa7fd429fa42c4eaf4ae5edc6a58d70e9d7f33))

- Removed `prepend` method attribute from `CifIO`, as it won't be needed with expanding...
  ([`ccf22e7`](https://github.com/Baharis/hikari/commit/ccf22e748bf777ff494aba20a958edfd6ae1cdf7))

- Removed `ustrip` function from `hikari.dataframes.cif` for the sake of including `uncertainties`
  later
  ([`c353a46`](https://github.com/Baharis/hikari/commit/c353a46b48b0b8400abbf5b9b9fc86adfc63b333))

- Removed dangerous setters from `BaseFrame` (after setting other pars wouldn't update)
  ([`0f1c929`](https://github.com/Baharis/hikari/commit/0f1c929d37e33181fb2b78e7541dc5303e24000f))

- Removed misleading if-name-main line in `hikari.dataframes.base` used for testing
  ([`c4a216e`](https://github.com/Baharis/hikari/commit/c4a216e14ef6a223b1ceafc2677c183f5540d9b2))

- Removed misleading if-name-self in `hikari.dataframes.cif`
  ([`55b6caf`](https://github.com/Baharis/hikari/commit/55b6caf7a1abf3a15649624780b3f5ae203dac47))

- Removed old TODO; CifFrame will exist besides PyCifRW
  ([`977739c`](https://github.com/Baharis/hikari/commit/977739cbd5ca6cd1d5e3fdeb6f8c1551137908ad))

- Removed older `compare_adp` script
  ([`c8eea1a`](https://github.com/Baharis/hikari/commit/c8eea1a615bf1291a4bbc0e07ed6be0695a06c2b))

- Removed project `requirements.txt`
  ([`f76d263`](https://github.com/Baharis/hikari/commit/f76d2639df6263d56645fc56b87440fc09e56d13))

- Removed trailing spaces from `.readthedocs.yaml`
  ([`19c781d`](https://github.com/Baharis/hikari/commit/19c781d317e8427d4007f45e4bc674c48f95c9d1))

- Removed whitespace padding from `BaseFrame` and `UBaseFrame` to increase maintainability
  ([`ec267d9`](https://github.com/Baharis/hikari/commit/ec267d90389e44b1f77133af8db869e15400e921))

- Renamed `compare_adps` script to more fitting `calculate_similarity_indices`
  ([`cba17e0`](https://github.com/Baharis/hikari/commit/cba17e0cd899d5ea6ae172f2cfc27351dcbc1855))

- Renamed ambiguous `DataBuffer.parse` to `DataBuffer.add_word`
  ([`f15c948`](https://github.com/Baharis/hikari/commit/f15c9487a94ddf5d477660f34b2331174c158b70))

- Script `calculate_similarity_index` now uses unit cell parameters from both cif blocks.
  ([`4777f30`](https://github.com/Baharis/hikari/commit/4777f30b115bbde6bab39348185d1c40ce7cc47a))

- Setup now correctly imports uncertainties package, making wheels and documentation - hopefully -
  work
  ([`e521f27`](https://github.com/Baharis/hikari/commit/e521f27997d8929f19e3e48a1ad0b17d43737ed3))

- Simplified `CifIO.DataBuffer` reset and added tab handling
  ([`50e8623`](https://github.com/Baharis/hikari/commit/50e86232ed210542bf4c92cde13c9e63c650a521))

- Split `revert_whitespace` and regex compiler out of `substitute_whitespace_in_quotes` in `CifIO`
  ([`27f3c3a`](https://github.com/Baharis/hikari/commit/27f3c3abdf6e251d6f87d98b86c3e641465eebd5))

- Testing sphinx integration of python typing for `euler_rodrigues_matrix` function
  ([`a0f0825`](https://github.com/Baharis/hikari/commit/a0f08254e49c8b00865c4264fc0f8744d4ba420a))

- Tests now use NaCl files in test folder instead of resources
  ([`2413ac7`](https://github.com/Baharis/hikari/commit/2413ac7f0ff6c187951cf023b2a309073deb620d))

- Turned off pdf/epub generation to prevent LaTeX issues
  ([`98a32d0`](https://github.com/Baharis/hikari/commit/98a32d025d021c2ac9342b59deb2fe1de940d4d7))

- Update codecov.yml
  ([`00673d0`](https://github.com/Baharis/hikari/commit/00673d0995bea2ca3730343986084eebe041ca19))

- Update codecov.yml
  ([`faa3785`](https://github.com/Baharis/hikari/commit/faa3785f6ddf1212a9b2834e003aea5a0a197620))

- Update codecov.yml
  ([`047755b`](https://github.com/Baharis/hikari/commit/047755b4cc8530b4abdf9050185ddcfb2caa5bcb))

- Update codecov.yml
  ([`9e138b4`](https://github.com/Baharis/hikari/commit/9e138b42763b768952abda1c06f267fadbeb583b))

- Update codecov.yml
  ([`26ff4a1`](https://github.com/Baharis/hikari/commit/26ff4a15cebf34bafca976880408122e69c6647e))

- Update README.md
  ([`4c762e1`](https://github.com/Baharis/hikari/commit/4c762e1ecb3d3e3c630abc0c6a0ecbb64aa08a12))

- Updated wrong link in `README.md`
  ([`86e398c`](https://github.com/Baharis/hikari/commit/86e398cef4a80dda638945ffbfa6f6d91174e5cd))


## v0.1.5 (2022-05-28)

### Other

- `angularheatmapartist.focus` is now a tuple instead of dict.
  ([`32f38c8`](https://github.com/Baharis/hikari/commit/32f38c8414bda2019680ad804d0bbbb391734782))

- `angularpropertyexplorer`s and consequently `r1_map` and `potency_map` now also correctly accept
  orientation vector
  ([`2f72a0c`](https://github.com/Baharis/hikari/commit/2f72a0c59adb65cc0d481ae53adc564da02d402c))

- `hikari.scrips.r1_map` and `hikari.scrips.potency_map` bug-fixed and splightly simplified
  ([`61e8f95`](https://github.com/Baharis/hikari/commit/61e8f954d02c7394587d043223f8728cf64752d9))

- `hikari.scrips.r1_map` moved to the scripts namespace and made to use angular heatmap artists
  ([`f663cef`](https://github.com/Baharis/hikari/commit/f663cefdcf0890d0907b630b8736cd64b3871e50))

- `matplotlibangularheatmapartist` now properly plots (100), (010), (001) arrows and focus points
  ([`30b3a58`](https://github.com/Baharis/hikari/commit/30b3a580d2c5eede596c45fa9865bfa0cf72c7d7))

- `potency_map` now draws focus points correctly in trigonal and hexagonal systems
  ([`9f28ad2`](https://github.com/Baharis/hikari/commit/9f28ad287a52316ac42b9773a8c820c6e22c4fc7))

- `potency_map` now logs average and quartiles of potency based on map data
  ([`342ceeb`](https://github.com/Baharis/hikari/commit/342ceeb268e8bd7f74fe685e33dece8de1061da1))

- `r1_map` and `potency_map` scrips were adapted to utilise `hikari.utility.Interval`
  ([`64509d0`](https://github.com/Baharis/hikari/commit/64509d0a2ce5d9d0ed8e932e9784b6eb99b58a36))

- `r1_map` now correctly determines trim area using radians, not degrees
  ([`09bdf32`](https://github.com/Baharis/hikari/commit/09bdf32c7494c0141f0149c39afd4cd3ffbd2eb0))

- Added `Interval` class docstring
  ([`127d11a`](https://github.com/Baharis/hikari/commit/127d11a610e2b32b5ddee43338dae7be4735d81f))

- All scripts use path created by `make_abspath` instead of raw strings
  ([`71c0bff`](https://github.com/Baharis/hikari/commit/71c0bff7aa558cdeee7f21668dd9756f8f2a6f47))

- Angularpotencyexplorer now correctly uses keys `potency` and `R1` instead of `cplt` and `r1`
  ([`767c460`](https://github.com/Baharis/hikari/commit/767c460fb1a045d58521b5e295e5284db8c39e81))

- Border lines are now drawn correctly in `gnuplot_angular_heatmap_template`
  ([`e1a2288`](https://github.com/Baharis/hikari/commit/e1a22888c6da626997d75db6bca9c8f681ef9178))

- Bumped hikari version to 0.1.5
  ([`f3dc57b`](https://github.com/Baharis/hikari/commit/f3dc57b2fc42dfb199779e2455f93131ba14fa3f))

- Changed `property_theoretical_minimum` and `..._maximum` to abstract in `AngularPropertyExplorer`
  ([`4df899f`](https://github.com/Baharis/hikari/commit/4df899fb9b384fa61f0ac9069b664f489f136f4c))

- Changed a number of attributes to properties in `AngularPropertyExplorer` (but still not
  functional)
  ([`890b886`](https://github.com/Baharis/hikari/commit/890b886f5d7368ae303b53079f642f28c81e75ba))

- Changed explore to abc method in `AngularPropertyExplorer`
  ([`3bc1cd0`](https://github.com/Baharis/hikari/commit/3bc1cd0d5a4248482df2e6ac0285176ff2ab9391))

- Created a `AngularPropertyExplorerFactory` for unified handling of angular property (r1, potency)
  maps
  ([`5eec073`](https://github.com/Baharis/hikari/commit/5eec0736ed006e8c050bb66e5921aa3f7cf2053a))

- Directional axes on `MatplotlibAngularHeatmapArtist` are now properly labelled (100), (010), and
  (001)
  ([`43e8329`](https://github.com/Baharis/hikari/commit/43e8329458444bdd2c298c3c86b2b52a112c4fe1))

- Drafted first version of `AngularPotencyExplorer` and `AngularR1Explorer` based on
  `AngularPropertyExplorer`
  ([`969b672`](https://github.com/Baharis/hikari/commit/969b6724ee8a07cd8a7d412c3daa8c9c548b3ce9))

- Fixed wrong `z_axis` setter in `AngularHeatmapArtist`
  ([`22bb896`](https://github.com/Baharis/hikari/commit/22bb8963dc638c8103f2370803a84e2d42ec7165))

- Further changes (property limit, class settings) introduces to still-wip `AngularPropertyExplorer`
  ([`cfe8393`](https://github.com/Baharis/hikari/commit/cfe8393203913a4e2369db8d500232cfb50acb76))

- Further minor docstring and naming changes of scripts in `potency_map.py`
  ([`067ab6c`](https://github.com/Baharis/hikari/commit/067ab6c5d1cec94e8f2e2f161d9502289448d393))

- Generated documentation (sphinx make html/clean) for hikari version 0.1.5
  ([`f0567c4`](https://github.com/Baharis/hikari/commit/f0567c40076501c72bc6106fa7ca9f7406148e8e))

- Groups and limits in `AngularPropertyExplorer` are now set by space_group setters, making code
  almost as fast as previously
  ([`2373929`](https://github.com/Baharis/hikari/commit/2373929131c7792cd463907585db20afeb50b24d))

- Hack-fixed matplotlib Axis3D aspect ratio by plotting invisible [1,1,1] and [-1,-1,-1] points
  ([`ff16db4`](https://github.com/Baharis/hikari/commit/ff16db46555b793df79bc37b19d843a1c7c5861c))

- Histogram made constant-width and map moved left in `gnuplot_angular_heatmap_template`
  ([`501fd5e`](https://github.com/Baharis/hikari/commit/501fd5edf314b1bb6eaedec54b55305178a3af00))

- Implemented `cartesian2spherical` and `spherical2cartesian` conversion functions
  ([`3f46d3c`](https://github.com/Baharis/hikari/commit/3f46d3c694e24edc872e5729dbaca8b88731edd2))

- Implemented `histogram` option in `potency_map`, for now does not draw
  ([`cabf6a2`](https://github.com/Baharis/hikari/commit/cabf6a21a3bfee232a31a9069823e39f80f3526c))

- Implemented a `weighted_quantile` in `utility.math_tools` which approximates quartiles in weighted
  distribution
  ([`e7dfdc4`](https://github.com/Baharis/hikari/commit/e7dfdc45ec5fad407d1e762dc520550e7a9663c8))

- Implemented method `arange` in `hikari.utility.Interval`
  ([`05d87e6`](https://github.com/Baharis/hikari/commit/05d87e65ca3a02f88c73b5f4740345175838952a))

- Implemented method `comb_with` in `hikari.utility.Interval`
  ([`a14206c`](https://github.com/Baharis/hikari/commit/a14206c6dbce10226a3c2f3bf7d5f7c07e798bcf))

- Implemented method `mesh_with` in `hikari.utility.Interval`
  ([`4779be3`](https://github.com/Baharis/hikari/commit/4779be385d21cfbc34941f5f24c35a101d8f53b3))

- Implemented unit tests for `hikari.utility.Interval`
  ([`8cacef0`](https://github.com/Baharis/hikari/commit/8cacef0373c85118cf55008d30e4bfc25ea231f3))

- In `MatplotlibAngularHeatmapArtist` made a bit more space for colorbox key
  ([`aff4037`](https://github.com/Baharis/hikari/commit/aff403717257ce0922c1e672233843fb604357b6))

- Included histogram-making in `gnuplot_angular_heatmap_template`, fixed pathing issues
  ([`742eb8b`](https://github.com/Baharis/hikari/commit/742eb8b0f3c826f2ed1d2c09d227ac79491b3c75))

- Introduced `line_subset` module from `pruby.utility.line_subset`
  ([`7c9517f`](https://github.com/Baharis/hikari/commit/7c9517fc42cd32c152522ba95ed40d73e09ff694))

- Merge pull request #15 from Baharis/paths
  ([`66f1450`](https://github.com/Baharis/hikari/commit/66f14509a7ac09c2a47934549436ae5b8ae6c97e))

All scripts use path created by `make_abspath` instead of raw strings

- Merge pull request #16 from Baharis/potency_map
  ([`5f64466`](https://github.com/Baharis/hikari/commit/5f6446604e881c5e7439a09295903f7e9f3e22d0))

Potency map

- Merge pull request #17 from Baharis/interval
  ([`3dc913f`](https://github.com/Baharis/hikari/commit/3dc913f2ba945ebabddb93ff9fd439b8e8de7106))

Interval

- Merge pull request #18 from Baharis/angular_explorer
  ([`a1b5785`](https://github.com/Baharis/hikari/commit/a1b57857ada0cf92962bfc1dabd3d52949909542))

Angular explorer

- Merged all `set_` methods` into one `set_up`, registered `AngularPropertyExplorers` in Factory,
  and replaced old `r1_map` and `potency_map` with new ones
  ([`a59f1de`](https://github.com/Baharis/hikari/commit/a59f1de52635bd9fc86ea8b1b4be9630d973f459))

- Minor clarity improvements after introducing `spherical2cartesian` function
  ([`575d5ee`](https://github.com/Baharis/hikari/commit/575d5ee1d8bee2d92a3b0aa3cd3988396c0209b7))

- Minor docstring and changes, unused items removes in `hikari.utility.artists` and
  `hikari.scripts.hkl_potency`
  ([`122ccf5`](https://github.com/Baharis/hikari/commit/122ccf5344f1434153f3bd9cb3e35ac2af088a45))

- Minor docstring and naming changes of scripts in `potency_map.py`
  ([`c79c0cf`](https://github.com/Baharis/hikari/commit/c79c0cfb924c126614f3a3cf451bf804c4579926))

- Modified `MatplotlibAngularHeatmapArtist` to plot larger maps (although histogram it too
  complicated to bother)
  ([`f249f4f`](https://github.com/Baharis/hikari/commit/f249f4ff61a1ba0fed3c7d7b7f1dff5ba5d029d7))

- Modified and merged previously imported `LineSubset` and `LineSegment`
  ([`109232e`](https://github.com/Baharis/hikari/commit/109232e00d1618dc56633c7c89ba228cfc2523c3))

- Modified recursive `_min` and `_max` functions used by `LineSegment`
  ([`20e5379`](https://github.com/Baharis/hikari/commit/20e53795cf201221f873ad253e5194e57551778e))

- Move `utility.artists` to `artist_factory` and unified plotting in `AngularPropertyExplorer`s
  ([`390cef4`](https://github.com/Baharis/hikari/commit/390cef4f2c64b5b39611e0b5befb641d65b503ef))

- Moved `PG` and `SG` declarations from `resources` to `symmetry` to avoid circular imports
  ([`fcce0c6`](https://github.com/Baharis/hikari/commit/fcce0c6bc40b27e023298c54dd26d5baf2fcd6f3))

- Moved gnuplot spherical heatmap plotting capability to `hikari.utility.artists`
  ([`9d1bfef`](https://github.com/Baharis/hikari/commit/9d1bfef4077eb02d78b4dc3e6b236f9eb3b3c16a))

- Moved matplotlib angular heatmap plotting capabilities to
  `hikari.utility.MatplotlibAngularHeatmapArtist`
  ([`f7a7d63`](https://github.com/Baharis/hikari/commit/f7a7d63dbec1f5da4e286790f9acf1d5226e035f))

- Potency map accepts argument `orientation` and appropriately plots `focus` point on the map
  ([`a900099`](https://github.com/Baharis/hikari/commit/a9000990bf9f2384dee7460fb5484649a01b4ba7))

- Removed old documentation
  ([`899d5a0`](https://github.com/Baharis/hikari/commit/899d5a0e94ad85a00f27b6b585635d7fb31ad7f5))

- Renamed `cartesian2spherical` and `spherical2cartesian` to `cart2sph` and `sph2cart`
  ([`d661340`](https://github.com/Baharis/hikari/commit/d661340fedc0f4d029a385e516e9fed4e08e95c4))

- Renamed `cplt_map_template` to `gnuplot_angular_heatmap_template`
  ([`b2824a7`](https://github.com/Baharis/hikari/commit/b2824a719b3629427be47a53e8785af735b8d926))

- Reorganised imports in `hikari.utility.__init__` to prevent import errors
  ([`00f59cb`](https://github.com/Baharis/hikari/commit/00f59cbec473626854f92b3e59a4b236801d8777))

- Simplifications and style changes in `gnuplot_angular_heatmap_template`
  ([`f9f0b63`](https://github.com/Baharis/hikari/commit/f9f0b6328f563cab8a956775b6a4136384bcccd4))

- Simplified `_min` and `_max` functions to be at most 1-deep recursive
  ([`657be28`](https://github.com/Baharis/hikari/commit/657be28ef9f2eb94091b6d9f11ed4017f6932d90))

- Simplified `potency_map` using new spherical2cartesian function
  ([`de94025`](https://github.com/Baharis/hikari/commit/de94025be51b3bbd9d4a6529ec3e4103adbf101b))

- Spherical heatmap now draws (100), (010), and (001) labels next to arrows
  ([`acdce23`](https://github.com/Baharis/hikari/commit/acdce23eacce082f066fd0894b028688ee5d1343))

- Started creating a `AngularPropertyExplorer` for handling of angular property (r1, potency) maps,
  not yet functional
  ([`168d0df`](https://github.com/Baharis/hikari/commit/168d0dfc4a7078049a8888a96e0b0f758763016e))

- Tested and fixed new potency_map2 and r1_map2 based on `AngularPropertyExplorer`
  ([`d88faf9`](https://github.com/Baharis/hikari/commit/d88faf9fefe4006838af07cca21840eb3ddaad2c))

- Triclinic and monoclinic systems have half and quadrant drawn.
  ([`e92db08`](https://github.com/Baharis/hikari/commit/e92db08a41533a705b74d2bdca290023fb131a78))

- Updated `__main__` script contents in `hkl_completeness.py` in `hikari.scripts`
  ([`3375bf5`](https://github.com/Baharis/hikari/commit/3375bf5f40a611824a11cf67814421bcd156f8ad))

- Updated `potency_map.py` and `r1_map.py` `__main__` scripts and added basic `r1_map` docstring
  ([`4cebc36`](https://github.com/Baharis/hikari/commit/4cebc36ddb7dc6f949094be490fa954e785086b2))

- Updated `potency_vs_dac_opening_angle` for potency and Python3.6 syntax
  ([`cec0efd`](https://github.com/Baharis/hikari/commit/cec0efd32627fdc7e4df83fa69186c05f4507a5c))

- Updted Interval unit tests to use modern `np.array_equal` instead of depreciated `(A==B).all()`
  syntax
  ([`569d22f`](https://github.com/Baharis/hikari/commit/569d22fa0a6a7953687def50f9d54d8763b10c2a))

- Upgraded and renamed `_assert_iterable_length` to `_assert_is_iterable` in `Artist`
  ([`d77332e`](https://github.com/Baharis/hikari/commit/d77332e813b9ea62aa42a5d239e575fef57c4a26))

- Vectorised `orientation_weight` to speed up `AngularPropertyExplorer`
  ([`71ce337`](https://github.com/Baharis/hikari/commit/71ce337c1e7f0e4b7469883c22cbb2b4b8c7a2b9))


## v0.1.4 (2022-04-20)

### Other

- Added and checked SFactors for Na+ and Cl- in NaCl, improvement only on low-res
  ([`37309f1`](https://github.com/Baharis/hikari/commit/37309f1ba30ba531d0e21ffd8abc3aa16de6f182))

- Added atomic_form_factor function to hikari.utility.chem_tools. Import errors still present.
  ([`e61a559`](https://github.com/Baharis/hikari/commit/e61a559718db5eaf624855410d1c633d83106b2e))

- Added correct orientation handling to the structure factor formula in ResFrame
  ([`c9e7bb6`](https://github.com/Baharis/hikari/commit/c9e7bb61a910264a0b25050c84f18fcd601c4982))

- Added docstrings to ResFrame and calculate_sample_form_factors script
  ([`06d1265`](https://github.com/Baharis/hikari/commit/06d126584cb2fe4e860d02ecc861bda055301be6))

- Added split_atom_label function to utility.chem_tools
  ([`f8742c8`](https://github.com/Baharis/hikari/commit/f8742c89ff4c690e354224ba212901d7f702d10a))

- Added temporary script to find best formula for temperature factor
  ([`0a7e7f1`](https://github.com/Baharis/hikari/commit/0a7e7f1c825f6af2353806448ba39bed7e9cebd0))

- Added temporary temp_calculate_other_form_factors function to fcf scripts
  ([`bbb40f3`](https://github.com/Baharis/hikari/commit/bbb40f3947d5892b5eed4d83d30f46a671968603))

- Added two rough methods to approximate structure factors, inaccurate for large sintl
  ([`1faa66e`](https://github.com/Baharis/hikari/commit/1faa66e0378b0daae71f6b8922f6ec4157cf1aba))

- Added Xray_atomic_form_factors to hikari.resources
  ([`402c601`](https://github.com/Baharis/hikari/commit/402c6012e36804dc8395e92fb97fc32af8f2631f))

- Changed formula to the one from Grosse-Kunstleve (2002)
  ([`6268fe9`](https://github.com/Baharis/hikari/commit/6268fe9de70cad04cee3272756c45006369145d4))

- Corrected Uij order as given in res file. Transformed Uijs are correct.
  ([`c35d1cb`](https://github.com/Baharis/hikari/commit/c35d1cbe16ac96eeaeced5bb88c65dcd9c66731a))

- Fixed the formula for atomic form factor, again
  ([`fc3d251`](https://github.com/Baharis/hikari/commit/fc3d251d7b208de0546b6b14ec3b6542023e1c74))

- Group.is_symmorphic is a bit slower, but works even if generators are unspecified or
  unconventional
  ([`0d4b4e3`](https://github.com/Baharis/hikari/commit/0d4b4e36d8c2437db471bf46374bd345f92a33b7))

- Group.name appends now +/- sense to axes. Inconsequent with ITC in cubic system, where the concept
  of axis direction is complicated
  ([`577fb1a`](https://github.com/Baharis/hikari/commit/577fb1a15a892a9a3a7f02d430c748629d7622bc))

- Hexagonal R-centred cells should are now called correctly
  ([`3c6e9ff`](https://github.com/Baharis/hikari/commit/3c6e9ff9671e8e94a6cedbe9d156d235b5306871))

- Introduced dummy LstFrame to look for value of R1
  ([`145ed9a`](https://github.com/Baharis/hikari/commit/145ed9a241078541dba92f8190e75c6609a1c3d5))

- Merge pull request #13 from Baharis/Rmap
  ([`287f0b5`](https://github.com/Baharis/hikari/commit/287f0b54d484b8701356071627ad235600f9bd33))

Rmap

- Merge pull request #14 from Baharis/res_frame
  ([`f0081ee`](https://github.com/Baharis/hikari/commit/f0081eea332e87b6682256efe95c796dcf72679b))

Res frame

- Morphed temps into a bare-bones script for fcf prediction
  ([`b738828`](https://github.com/Baharis/hikari/commit/b7388288ddef6026a4b2f0b7b07a1aca6482aecb))

- Prepared a simple script to map R1 using shelxl
  ([`904d8da`](https://github.com/Baharis/hikari/commit/904d8da68f0ad5ee99e7efa7dc39bd39fbb68eb8))

- Prepared setup.py file and introduced as hikari-toolkit to pypi
  ([`1043ef2`](https://github.com/Baharis/hikari/commit/1043ef2a08f38adb5381d3ef640fd2185ac31766))

- Slightly rewritten existing read method of ResFrame to understand syntax
  ([`1826a24`](https://github.com/Baharis/hikari/commit/1826a24a73659278eac03b55c96bf6e837a23337))

- Source of Xray_atomic_form_factors added to csv. Some circular import errors found.
  ([`22842e3`](https://github.com/Baharis/hikari/commit/22842e3ce1de9452af298c34f4cfa912d25d8527))

- Split atomic_form_factor and temperature_factor into separate methods
  ([`dbd13a3`](https://github.com/Baharis/hikari/commit/dbd13a3345b7952c16a16b70345ab46d0fe7951e))

- Symmop.orientation now utilises a direct formula instead of invariants. Rotations still have funny
  names in cubic system.
  ([`e6a9db8`](https://github.com/Baharis/hikari/commit/e6a9db83f62e9b68169cedd5e40aef649e37fcea))

- Update README.md
  ([`eb98381`](https://github.com/Baharis/hikari/commit/eb9838177c4119c6b3889fca3140353872b7a7ba))

- Vectorised atomic_form_factor, temperature_factor, and form_factor methods
  ([`a4da118`](https://github.com/Baharis/hikari/commit/a4da11878df2eac2176fd3d6beafba803b98006a))


## v0.1.3 (2022-01-22)

### Other

- Added additional simple tests for SymmOp
  ([`de7ff2d`](https://github.com/Baharis/hikari/commit/de7ff2d7149df59dc7d02c83afa5bcc6d743fc3c))

- Added comparison methods to hikari.symmetry.Group
  ([`a01297e`](https://github.com/Baharis/hikari/commit/a01297e8c0d4d264ce88a4da2e858eb0ba52910a))

- Added crystallographic matrices A and G to BaseFrame
  ([`79eefd4`](https://github.com/Baharis/hikari/commit/79eefd4f04f6ce7e4c75019074347717a05975dd))

- Added hkl.msd style file to resources, printed with each res
  ([`e1e4c0b`](https://github.com/Baharis/hikari/commit/e1e4c0b60735ab138790d78fc67897a0b3dcb45d))

- Added human&computer-readable representation of SymmOp
  ([`7660c47`](https://github.com/Baharis/hikari/commit/7660c47899c15915f64df83da358171e58d320c2))

- Added sample unittest for hikari.symmetry.operations
  ([`3d6ea13`](https://github.com/Baharis/hikari/commit/3d6ea1387db80e7083538ca7fcdaa750cea6400b))

- Added simple tests for hikari.utility
  ([`9dc37bb`](https://github.com/Baharis/hikari/commit/9dc37bb9c8d9f6d888e30aca38e59f65caf64de7))

- Based DAC definition on matrix A
  ([`340264f`](https://github.com/Baharis/hikari/commit/340264f9bb3dd61fc0a63388f919fc80ae7a2b09))

- Changed definition of BaseFrame.v to determinant of matrix BaseFrame.A
  ([`597f3ec`](https://github.com/Baharis/hikari/commit/597f3ec2249a4336fa653294195476cb59fa647e))

- Characteristic radiation wavelengths moved from HklFrame.la to json in resources
  ([`30ed90c`](https://github.com/Baharis/hikari/commit/30ed90c7112da201335ce06c38da22f76994ffd3))

- Cplt_map_template.gnu renamed and moved to new directory resources
  ([`f4e068e`](https://github.com/Baharis/hikari/commit/f4e068e10c182a55b2d362a79963cf0a55b43878))

- Created fully vector dac_count unfortunately has enough memory only for small system, so reverting
  to hybrid approach
  ([`d07a5b0`](https://github.com/Baharis/hikari/commit/d07a5b0470db506e5f277c68c8d3792da74046c4))

- Dac_count of HklFrame accepts only one vector now, but is 2x faster then previous solution due to
  vectorisation
  ([`b168530`](https://github.com/Baharis/hikari/commit/b1685304c022561b6f8db8bf9e7e93e2ea317db0))

- Dac_potency_around_axis script changed to utilise new rotation from math_tools
  ([`91f9ee7`](https://github.com/Baharis/hikari/commit/91f9ee7cb6b6bf464f9283c9dfca0d816da7998f))

- Dacs_count now splits input into chunks not to exceed hikari.MEMORY_SIZE
  ([`3f39710`](https://github.com/Baharis/hikari/commit/3f39710a8d2e92cdbd2b7f3ba226a152f870a3bd))

- Defined rotation matrix in math_tools using Euler-Rodrigues formula
  ([`29f8d23`](https://github.com/Baharis/hikari/commit/29f8d238b17795ed3f77235526aea7199e860452))

- Duplicated potency_map function for further edition
  ([`c9d6849`](https://github.com/Baharis/hikari/commit/c9d6849f15bf9211fa31126670cc9c7b500ede69))

- First draft of fully vectorised dac_count function for maps
  ([`b528fdb`](https://github.com/Baharis/hikari/commit/b528fdb1f1931f95083adc8b72fc05024e5f001c))

- Fixed (temporarly?) circular import issue in hikari.dataframes.hkl
  ([`246715e`](https://github.com/Baharis/hikari/commit/246715e02359dc17c913d7a14308254ba62e881d))

- Fixed extinct method of HklFrame broken after last name changes
  ([`cf37198`](https://github.com/Baharis/hikari/commit/cf37198e70ff4988cfe379bb264997323e86c468))

- Fixed HklFrame.HklIo fixed parser, which for some bizarre reason returned None
  ([`1bb4d3b`](https://github.com/Baharis/hikari/commit/1bb4d3bc341632fa49e3ea27b44bda39befd3e96))

- Fixed x/y/z definitions in new converter and changed it to default artist in HklFrame
  ([`2c57a11`](https://github.com/Baharis/hikari/commit/2c57a11be845a9d1ba8f6abe118b838d4e0a0ad6))

- Hikari.dataframes.hklframes().extinct() now correctly uses Group().operations
  ([`7c688ab`](https://github.com/Baharis/hikari/commit/7c688ab73cb1b4e4e15df25fe77dd9d9abe9bfc3))

- Hkl format dictionaries moved to resources and supplied directly as dict
  ([`1bf71e5`](https://github.com/Baharis/hikari/commit/1bf71e54604df0bb23af8182af20ddb73a3e0adb))

- Hklframe._in_dac further optimised by calculating norm directly
  ([`29a8931`](https://github.com/Baharis/hikari/commit/29a893107ad53c9e405b543c430e102eb9bf4903))

- Hklframe._in_dac removed, all references rerouted to vectorised _in_dacs
  ([`1375add`](https://github.com/Baharis/hikari/commit/1375add8df69e1c3b4dd3df4df2c1421d05abf2d))

- Hklframe.fill simplified to use reshape instead of nested concatenates
  ([`a77ee4e`](https://github.com/Baharis/hikari/commit/a77ee4eb8e3808e3a1422310da3a67a165c4bfd7))

- Hklframe.merge docstring abbreviated
  ([`c4aacd9`](https://github.com/Baharis/hikari/commit/c4aacd9b1e26c5554adfd632840a8266a6539798))

- Hklframe.stats now returns correctly cplt and I/si
  ([`e94ffab`](https://github.com/Baharis/hikari/commit/e94ffab97a5f61af30eda091af9070a97c09db68))

- Implemented 1st draft of HklToResConverter to replace HklArtist
  ([`61d3112`](https://github.com/Baharis/hikari/commit/61d31126245be2397fd1f38fc72681f895872d67))

- Implemented proper str, repr and eq methods for Group
  ([`ba52cb6`](https://github.com/Baharis/hikari/commit/ba52cb6a205b07b6819227099eee14e90403ff01))

- Introduced more advanced tests for HklFrame, for trim, place etc
  ([`a87828b`](https://github.com/Baharis/hikari/commit/a87828b00305f48fe72c6213ac5811dbbcb92a0d))

- Introduced simple tests for all utility functions and symmetry.operations
  ([`e5f72f9`](https://github.com/Baharis/hikari/commit/e5f72f9ffc4ef6dc7de32485d993dbeb5653b85b))

- Introduced simple tests for HklFrame
  ([`d43b6e5`](https://github.com/Baharis/hikari/commit/d43b6e5a6b8a32fde7d6758e84ac4271fa655981))

- Introduced simple tests to control BaseFrame algebra
  ([`4366031`](https://github.com/Baharis/hikari/commit/4366031deecafe2aadd36041f1f51499e369a910))

- Merge branch 'master' into unittest
  ([`2dd837c`](https://github.com/Baharis/hikari/commit/2dd837c0db385c5c5a5eba0931f1538523cdd386))

- Merge pull request #10 from Baharis/fibonacci
  ([`53647f6`](https://github.com/Baharis/hikari/commit/53647f6f00c2a0481352ffd59fb04c4140e3bc5f))

Fibonacci sphere algorithm

- Merge pull request #11 from Baharis/scripts_simplified
  ([`ad552ad`](https://github.com/Baharis/hikari/commit/ad552ad7f3b9796eefad1c94016f4acf6311a1c3))

dacs_count debugging

- Merge pull request #12 from Baharis/frame_tests
  ([`22c812f`](https://github.com/Baharis/hikari/commit/22c812f98347744f4c8aaf00cb5c2f68ce4019d4))

Frame tests

- Merge pull request #3 from Baharis/resources
  ([`51e0f1b`](https://github.com/Baharis/hikari/commit/51e0f1b8536d31e0a5caa77e4c961e3eace9a475))

Resources management optimisation

- Merge pull request #4 from Baharis/representations
  ([`f7229d5`](https://github.com/Baharis/hikari/commit/f7229d50628ab4c7dc31b8afab49bd34b0204ffe))

Representations

- Merge pull request #5 from Baharis/hkl_artist
  ([`9f7ebc9`](https://github.com/Baharis/hikari/commit/9f7ebc910184db4a5e3fd35762254c4a341f20f6))

Replaced HklArtist with simpler HklToResConverter

- Merge pull request #6 from Baharis/stats
  ([`8742bdf`](https://github.com/Baharis/hikari/commit/8742bdf36dd3a23ef656d6dc8c3dd583bbea4c04))

Method `stats` of `HklFrame` used to return wrong completeness statistics for merged datasets. After
  this fix, `stats` not only calculates completeness correctly, but also "I/si(I)". This might fail
  whenever for some reason the file does not feature `I` or `si`.

- Merge pull request #7 from Baharis/hklframe
  ([`f9c40fb`](https://github.com/Baharis/hikari/commit/f9c40fbfc2367bc3035f8bd447e8c613057d2db0))

`hikari.dataframes.hkl.HklFrame` methods and doc-strings simplified, circular import fixed.

- Merge pull request #8 from Baharis/dac_count
  ([`ac79fc6`](https://github.com/Baharis/hikari/commit/ac79fc6af9d728060d28a5e6ec4f2757949f4f60))

dac method optimised and split into dac_count and dac_trim

- Merge pull request #9 from Baharis/unittest
  ([`9bcdb0a`](https://github.com/Baharis/hikari/commit/9bcdb0acfc2f996d0e29ece6e8d74ab9f412ba5e))

Simple unittest-s for utility functions and symmetry.operations

- Minor import optimisations
  ([`1617786`](https://github.com/Baharis/hikari/commit/16177864d887f8435c3cc3fb8e53cfc3ac79ae61))

- Minor modifications, potency_map 1 and 2 merged
  ([`0b9f36b`](https://github.com/Baharis/hikari/commit/0b9f36ba858cbbb1c09fde955705a4b89619b2c8))

- Moved priority rules from Group().auto_generated_rules to Group attributes
  ([`aacfb8b`](https://github.com/Baharis/hikari/commit/aacfb8b981186417e0347c0a79aee248cff99639))

- Point and space groups are now loaded in resources and only aliased in symmetry's init
  ([`282469c`](https://github.com/Baharis/hikari/commit/282469ceea0ec5b466ccdd7a4a7b7e47e9a8c66f))

- Point_groups.json and space_groups.json moved to resources
  ([`5a78d9a`](https://github.com/Baharis/hikari/commit/5a78d9af2c261468f3fa7bf19552c49fa41c0144))

- Point_groups.pickle and space_groups.pickle moved to resources
  ([`53dc41d`](https://github.com/Baharis/hikari/commit/53dc41d6ee749c5d6ee21cb9eb6446c38f1b46b9))

- Potency_map code clarity fixes
  ([`2d5bc8a`](https://github.com/Baharis/hikari/commit/2d5bc8afc962e33f3c0b24551eecbe02b5228693))

- Potency_violin_plot script changed to utilise new dacs_count
  ([`3f2e2aa`](https://github.com/Baharis/hikari/commit/3f2e2aa4e77761a0f9633426e57de7f77a6715f5))

- Potency_vs_oa script changed to be more precise and use new dacs_count
  ([`cef7470`](https://github.com/Baharis/hikari/commit/cef7470e5060b42ee5d45871ed6de778993032ba))

- Re-implemented fully-vectorised dacs_count added and removed in a previous commit
  ([`06dbcba`](https://github.com/Baharis/hikari/commit/06dbcba5110fd1a43ab6c603862f148de4ffecda))

- Re-implemented fully-vectorised potency_map script, which is 10-20% faster
  ([`0fa2266`](https://github.com/Baharis/hikari/commit/0fa226657f4bd2db9f7c15886dcdec7c98c010c6))

- Refactored HklFrame.duplicate to HklFrame.copy
  ([`150fb3b`](https://github.com/Baharis/hikari/commit/150fb3b6af2bdc8a5c38358ca8abb31e397441d2))

- Removed dac_count as redundant, simplified docstrings
  ([`d72ad28`](https://github.com/Baharis/hikari/commit/d72ad28acb9e3a47e7b0e3367c2ec1d616334bd0))

- Removed memory logging and other comments from dacs_count
  ([`f18569f`](https://github.com/Baharis/hikari/commit/f18569fe35ca06eca6a39fc1dbb311a81fc6c3f3))

- Removed old, more complicated hkl-to-res converter.
  ([`1b45f42`](https://github.com/Baharis/hikari/commit/1b45f42768223b58bda57945efc61efd494cef5a))

- Renamed script files, moved all potency methods to separate file
  ([`a71af82`](https://github.com/Baharis/hikari/commit/a71af822439058dda070df7899a4e3a698c89435))

- Save pickle moved to hikari.resources
  ([`8611057`](https://github.com/Baharis/hikari/commit/861105781617a5b457113e8ba50dad9b2119bf7d))

- Scripts modified to accept space group string or int instead of Group instance
  ([`1784514`](https://github.com/Baharis/hikari/commit/1784514a1671f725026c481b9a55e3d4193f191f))

- Shortened some docstrings, but still no bugfix
  ([`7159e75`](https://github.com/Baharis/hikari/commit/7159e753606683888b8dd34148516ba06288dfb9))

- Simplified and fixed previously-bugged BaseFrame._refresh_cell
  ([`ca0753e`](https://github.com/Baharis/hikari/commit/ca0753e969e727e8891648fd3f33ea1ab0c0746c))

- Simplified BaseFrame documentation, removed x/y/z attributes to introduce orthogonal vectors later
  ([`db31cd7`](https://github.com/Baharis/hikari/commit/db31cd798df1457cc9a660820fccfa7a8148c850))

- Simplified docstring of hikari.scripts.potency_map
  ([`84f4782`](https://github.com/Baharis/hikari/commit/84f47821b1fc942b1a9e6eb91a489eabef149891))

- Simplified HklFrame.stats docstring
  ([`90be715`](https://github.com/Baharis/hikari/commit/90be715fa2299bcf5dfcee152f84c51175d3afbe))

- Simplified some tests to control BaseFrame algebra
  ([`6a6256f`](https://github.com/Baharis/hikari/commit/6a6256f67e3f456d52150160a53234bdee708166))

- Single vector HklFrame.dac_count restored not to use caching, which was buggy
  ([`3b2d41f`](https://github.com/Baharis/hikari/commit/3b2d41f50cb49e90c4ac32d066ee255ec9871c22))

- Split HklFrame.dac into dac_count and dac_trim based on vectorised _in_dac, adapted potency_map
  ([`6914e25`](https://github.com/Baharis/hikari/commit/6914e250c7deef43ec24d0f17822e1576de01f17))

- Symmop's hash now uses hash of SymmOp's repr
  ([`3ce639c`](https://github.com/Baharis/hikari/commit/3ce639c032e4a68c5e0c824b416ec59747944ac6))

- Tested and found fully-vectorised method to be just as fast as bugged-heuristic looping dac_trim
  ([`0739c55`](https://github.com/Baharis/hikari/commit/0739c55bde95511595286db5872d765f0cecf3c3))

- Uniformity tests of fibonacci_sphere now uses larger sphere of 99k elements
  ([`4b4f5fd`](https://github.com/Baharis/hikari/commit/4b4f5fdf154d517d5f220e9e0a40f8b7b9f7fa1d))

- Utilised A_r in HklFrame.place()
  ([`df753b0`](https://github.com/Baharis/hikari/commit/df753b089d341bfb890f2815f2b05c8fdff970fa))

- Vectorised hikari.utility.math_tools.fibonacci_sphere using numpy
  ([`ff6c969`](https://github.com/Baharis/hikari/commit/ff6c96940bf2504f247f390f4fd42c4bcff60638))


## v0.1.2 (2021-11-04)

### Other

- Added and modified some properties of hikari.symmetry.group
  ([`40cd5b6`](https://github.com/Baharis/hikari/commit/40cd5b654a61a17c13c6f61faa853d039fa1318e))

- Added approximate group naming capability. Group.system now has property directions.
  ([`fa44d89`](https://github.com/Baharis/hikari/commit/fa44d899037f469f73dd82e6aa29bbf84dd77181))

- Added number to group (default 0). Re-pickled PG & SG dicts.
  ([`b5dd93c`](https://github.com/Baharis/hikari/commit/b5dd93c57be010fc35dec1abcbe8163319a603d8))

- Code minor corrections and simplifications in scripts.hkl
  ([`a8bc84c`](https://github.com/Baharis/hikari/commit/a8bc84c3fd17354465b300bbc4a5e8bd4ff6ce63))

- Fixed bugs: 3>-3 and occasional 21>2 rules in Group.auto_generated_name
  ([`18e3558`](https://github.com/Baharis/hikari/commit/18e35583396c7b1167303261176078c5878aa8f5))

- Group.name hidden as Group._generated_name as it might be useful later
  ([`eb5ddd2`](https://github.com/Baharis/hikari/commit/eb5ddd23ad597208f5d397c370494c24901c82f8))

- Hikari.utility.home_directory deleted, as new make_abspath makes it obsolete
  ([`890083d`](https://github.com/Baharis/hikari/commit/890083d2d710c6079f52f901c93190156678c678))

- Implemented Group.transform() to enable using undefined settings.
  ([`e7702cf`](https://github.com/Baharis/hikari/commit/e7702cf80ef0aac70da6864cd78895a362b29aee))

- Introduced/reintroduced point/space group json files.
  ([`8150895`](https://github.com/Baharis/hikari/commit/81508959c64272de91a210c3dcd46ed92a178586))

- Make_abspath now correctly returns string instead of Path object
  ([`dac3346`](https://github.com/Baharis/hikari/commit/dac3346714fb58998df4227ca10939590c9b4d6e))

- Make_abspath renamed, now uses pathlib and correctly expands '~' & '..'
  ([`ca1866f`](https://github.com/Baharis/hikari/commit/ca1866fee9ebd68ff6d03ef47dd259ec91bf433d))

- Merge remote-tracking branch 'origin/master'
  ([`126ac6a`](https://github.com/Baharis/hikari/commit/126ac6a68ddec5bd71a7a1a2a11ff5a32baef594))

- Output of pandas.describe() does not collapse horizontally anymore
  ([`2aae422`](https://github.com/Baharis/hikari/commit/2aae4225c033906cbc9a8091672e02e7c748a5c2))

- Output of pandas.describe() is now appended to potency_violin_plot() log file
  ([`619702d`](https://github.com/Baharis/hikari/commit/619702d15dc9efa293ecaf2ed97ea75c50ab5178))

- Palettes.py in utility simplified, unused sns palettes removed
  ([`524c128`](https://github.com/Baharis/hikari/commit/524c12840b7259fce46ba31622459aa44fe0b5e8))

- Removed 'uncertainties' from requirements.txt
  ([`0b9af79`](https://github.com/Baharis/hikari/commit/0b9af7903957c69b99ec4175b1f8b499b61a7534))

- Removed unused imports from CifFrame and ResFrame
  ([`24db35a`](https://github.com/Baharis/hikari/commit/24db35ac7806d4383dbca8ba8df2c6054b17e329))

- Simplified definitions and docstrings in hikari.symmetry.group
  ([`980d8c3`](https://github.com/Baharis/hikari/commit/980d8c3d671e629975b74996537440be10a93614))

- Typo fix in README.md
  ([`5d9b18d`](https://github.com/Baharis/hikari/commit/5d9b18db7f152eb0964aa92339d2ead9e864e8a3))


## v0.1.1 (2021-10-20)

### Other

- Abbreviated Group.CrystalSystem to Group.System Fixed wrong definition in Group.system() which
  made mmm monoclinic
  ([`b3ba4c9`](https://github.com/Baharis/hikari/commit/b3ba4c9a5212596f7507a5f0a5f931d69ddaaa76))

- Added legacy folder to .gitignore
  ([`7f3cc87`](https://github.com/Baharis/hikari/commit/7f3cc87ce6cf1b6b7134b0358914fa20940e9f73))

- Changed HklKeys.__equiv type from tuple to int, resulting in 1.5-3x speedup.
  ([`8ec4e46`](https://github.com/Baharis/hikari/commit/8ec4e46d7ad063c501f63f178aa89457db2326bf))

- In scripts.hkl renamed completeness_map to potency_map.
  ([`dc53830`](https://github.com/Baharis/hikari/commit/dc538300617e2d8989ca0d25462447d3af745c11))

- Minor code clarity corrections.
  ([`3b04e52`](https://github.com/Baharis/hikari/commit/3b04e521157833a071d3c5d670c51b869164a217))

- Moved fcf and potency_map scripts to separate files, heavily optimised and renamed
  completeness_violin_plot to potency_violin_plot, merged dac_point_group_statistics into
  potency_violin_plot
  ([`e2af144`](https://github.com/Baharis/hikari/commit/e2af144727f3959039f749b954601cc9cc71a827))

- Point groups are now read from point_groups.pickle
  ([`28cc3c9`](https://github.com/Baharis/hikari/commit/28cc3c9468210bdadf2a797efedaad7482df60f0))

- Removed legacy is#n functions from hikari.utility.
  ([`f1c10cf`](https://github.com/Baharis/hikari/commit/f1c10cff34b5ef451951b9ea19b5fb297be2c84f))

- Removed legacy subfolder from repository
  ([`0f00b54`](https://github.com/Baharis/hikari/commit/0f00b54ee251cbd4bb2a8d19c459ddfc363117ed))

- Removed trimming in Hkl.dac(), increasing speed when called many times.
  ([`7047e98`](https://github.com/Baharis/hikari/commit/7047e98b5fee347c1165a2a0b0ce82c1b23442d3))

- Space groups are now read from space_groups.pickle
  ([`4723c51`](https://github.com/Baharis/hikari/commit/4723c5160a738dae1b9d9324c1ac2dca088a7bfa))


## v0.1.0 (2021-10-13)

### Other

- A lot of script modifications, fixed wrong transformation in hkl dataframe, slightly modified map
  palettes
  ([`8499ded`](https://github.com/Baharis/hikari/commit/8499ded8e02e986f0d78b16dd8d15213a1e570fb))

- A lot of stuff I don't remember
  ([`c3d4b30`](https://github.com/Baharis/hikari/commit/c3d4b30abb0653d8d72857cd4c84ee2c30961ff9))

- Added crystal system; about to modify group/space group
  ([`3915e47`](https://github.com/Baharis/hikari/commit/3915e477a6713b344202caa536f8ebf9f663c7ef))

- Added documentation and increased overall clarity of symmetry module.
  ([`33c6d1e`](https://github.com/Baharis/hikari/commit/33c6d1e40d0ba5f9395fd16162cf0ded8dedff27))

- Added documentation and increased overall clarity of utility module.
  ([`d139c11`](https://github.com/Baharis/hikari/commit/d139c111c075df651678cb681b5496a89f2f1aa6))

- Added example folder, generalised rescale method of HklFrame
  ([`5fffc76`](https://github.com/Baharis/hikari/commit/5fffc76af86b64081931ea03fb52ec7befac86a6))

- Added extinction functionality
  ([`9d02c6f`](https://github.com/Baharis/hikari/commit/9d02c6fd586217fe5a2b4b4a1c6afc6bd170f054))

- Added free format read/write (untested)
  ([`3bc335d`](https://github.com/Baharis/hikari/commit/3bc335d2724617bc04411314cf4821500f7af6a5))

- Added new concept of symmetry operation management
  ([`b2b7561`](https://github.com/Baharis/hikari/commit/b2b756142adf4be481c8374854ac643bd9846c19))

- Added rotoinversion type and into methods
  ([`0f2a486`](https://github.com/Baharis/hikari/commit/0f2a486d35fec6d355ccb27bb7c96182b1e2416a))

- Added some sphinx files
  ([`cb4808c`](https://github.com/Baharis/hikari/commit/cb4808c3a8b7538fb5a572da3853ece2c49c6c4d))

- All SGs have been defined in json file
  ([`74ff2cc`](https://github.com/Baharis/hikari/commit/74ff2cc214878acad2a59218dfe45a17743d4c67))

- Alternative SG origins have been defined in json file
  ([`fe69c1c`](https://github.com/Baharis/hikari/commit/fe69c1c911e34d0aad7c2c4d939372226ee75733))

- Changed names of point groups in documentation For various scripts: > Fixed completeness to
  calculate unique reflections instead of all resymmetrified (cplt1 instead of 2) > Added
  legacy_cplt to calculate cplt1
  ([`f3c64b4`](https://github.com/Baharis/hikari/commit/f3c64b466ec82483ca11738522e8dbaef7944168))

- Changed paths in gnuplot to absolute Fixed transform to act as resymmetrify before Changed
  disc-transforming ops to chiral_ops
  ([`d15cf20`](https://github.com/Baharis/hikari/commit/d15cf2000cca0c8f14622aa44f0f95225b52c68d))

- Changed point group to dictionary Simplified symmetry structure Simplified wavelength description
  Extended 'edit_wavelength' method
  ([`36dab53`](https://github.com/Baharis/hikari/commit/36dab5308c418deb50f543551d448a5f0e825a95))

- Changed SGs to rely on ITC-A generator definitions SGs are now stored in json file
  ([`1f64765`](https://github.com/Baharis/hikari/commit/1f64765061d2676ab79f7c437584c993ce17ec1d))

- Completeness vs DAC orientation script added
  ([`718db7c`](https://github.com/Baharis/hikari/commit/718db7cb890f7d1344730a88d14bfc339d18c53b))

- Cplt map script for plotting in gnuplot simplified
  ([`f7f8190`](https://github.com/Baharis/hikari/commit/f7f8190efcc673f3e1aa93c67d72f11af37e0357))

- Created taskmasters, Added normalized direction vectors x,y,z/v,w
  ([`77a0ab8`](https://github.com/Baharis/hikari/commit/77a0ab886a88723f355e5ce9d0ee890240dc7e45))

- Defined nnp, simple statistical tools and axis-specific colormaps
  ([`388fe3c`](https://github.com/Baharis/hikari/commit/388fe3cd95bcc9f20ad00d06f576b25c75d5e5cb))

- Defined space groups 1-46
  ([`7f05615`](https://github.com/Baharis/hikari/commit/7f0561586d0a6fe17feae0505967482af5e3702f))

- Delete 1%2E1686607.pdf
  ([`eca521d`](https://github.com/Baharis/hikari/commit/eca521dd65c86c80ccd50ce8f7876d82f2e34ae6))

- Deleted scripts dir
  ([`9d8a899`](https://github.com/Baharis/hikari/commit/9d8a899b136859b674076f51e71a849d3de814d4))

- Diamond anvil cell completeness map drawing script introduced
  ([`916aaf0`](https://github.com/Baharis/hikari/commit/916aaf07afd5fc6986cd50c3881686d4681129cc))

- During DAC tests
  ([`a083b38`](https://github.com/Baharis/hikari/commit/a083b381503ead212d4320559544ae74e0c5b24e))

- Finished clearing files, removed old components, added manual file
  ([`862d4be`](https://github.com/Baharis/hikari/commit/862d4be28eea702ccb34e85250e77463d9472496))

- Fixed colormap in matplotlib to show actual colors Added "fix_scale" option to completeness map
  script/taskmaster
  ([`5e795cd`](https://github.com/Baharis/hikari/commit/5e795cddc1e1df541891ce1453433975cb54fbea))

- Fixed reciprocal definition
  ([`caedf1d`](https://github.com/Baharis/hikari/commit/caedf1d5eb417aa505cf8efd0c7522a83988f623))

- Fixed SettingWithCopyWarning and some legacy code from numpy
  ([`0c69296`](https://github.com/Baharis/hikari/commit/0c69296ceb252a9a42ae69222a27a6fd39dc58cb))

- Fixed, cleaned and documented new SymmOp
  ([`241981c`](https://github.com/Baharis/hikari/commit/241981c7c63e003ea5a62026756d32c3b3f40797))

- Generalized extinct to accept tuple or string with ':' Added separate 'find_equivalents' method
  ([`bdcf053`](https://github.com/Baharis/hikari/commit/bdcf05366b8ea841eb7a34b09e73aff5d7111d02))

- Hikari and pRuby imported
  ([`7b6c0f8`](https://github.com/Baharis/hikari/commit/7b6c0f80c172ea0eb0b42fb6521b80aaa469f865))

- Hkl calculating moved
  ([`57c5e3f`](https://github.com/Baharis/hikari/commit/57c5e3fc2448e742c9b613005139f77e832c923f))

- Hklio introduced to handle reading and writing
  ([`75acd2e`](https://github.com/Baharis/hikari/commit/75acd2e7948cf2e3b3b86c6e9479e6131796aa75))

- Hklres visualiser added
  ([`b5e3954`](https://github.com/Baharis/hikari/commit/b5e3954c1206d078617bf66fdbcd9edbda7b89f1))

- Hklres visualiser added
  ([`f05d6fa`](https://github.com/Baharis/hikari/commit/f05d6fa1e244d3cfeed87d834b2171db81f07ea9))

- I made this
  ([`7d366e1`](https://github.com/Baharis/hikari/commit/7d366e116eb2cb7ac1e2ec873730b50203cbf2c8))

- Initial commit
  ([`ed432e3`](https://github.com/Baharis/hikari/commit/ed432e3b58ea058b175c93717392bdf95005e450))

- License and readme added
  ([`c353c2d`](https://github.com/Baharis/hikari/commit/c353c2db8b4bd19fc802abc678ccf8663633e78c))

- License and readme added
  ([`4dcbefd`](https://github.com/Baharis/hikari/commit/4dcbefd93cf37f2b7e0d1a56b64db08cdc2d09ec))

- License and readme added
  ([`1932d6f`](https://github.com/Baharis/hikari/commit/1932d6fb97d964a903728a6ca724d0bd73d51123))

- License file updated for year 2021.
  ([`512822d`](https://github.com/Baharis/hikari/commit/512822d09cd8b5c85f2522189dc70ac9ca7a08d3))

- Made hkl frame accept fcf files Prepared script for BayCoNs and partially npps
  ([`d2195b6`](https://github.com/Baharis/hikari/commit/d2195b63384d34bf387d2c08099fa3771017db29))

- Merge pull request #1 from Baharis/dataframes
  ([`f84fc81`](https://github.com/Baharis/hikari/commit/f84fc819a75402263d206c491fb9aff539c87beb))

Dataframes

- Merge pull request #2 from Baharis/dataframes
  ([`0b879eb`](https://github.com/Baharis/hikari/commit/0b879ebf0fcb3b480dc908321c38bbd3a1c70af8))

Dataframes

- Merged 'resymmetrify' and 'transform' Fixed 'make_stats' for given symmetry
  ([`780fee1`](https://github.com/Baharis/hikari/commit/780fee11f5d0bd4e6a41f9c251fac059935989e8))

- Minor changes to new symmetry operation algebra
  ([`46c55b4`](https://github.com/Baharis/hikari/commit/46c55b40e6cf55c760a1811fe2e2f5eb18598d1d))

- Minor correction to check git connection
  ([`07d0678`](https://github.com/Baharis/hikari/commit/07d06785d249735efe08cb8843eccd5cfb4f4d3a))

- Minor naming corrections from kesshou to hikari.
  ([`6f2cfe0`](https://github.com/Baharis/hikari/commit/6f2cfe02ff3a67c1a3ccd3b0f35c401b653ea4d7))

- Move unused code to external folder with stubs, fix documentation errors
  ([`8d1973c`](https://github.com/Baharis/hikari/commit/8d1973c7bdc65bc0f89bf3f23fd08a4c5d291eea))

- Moved all group and symmetry functionality to new class
  ([`219d1bd`](https://github.com/Baharis/hikari/commit/219d1bde578d495b407efd31a2418bbecda699e8))

- Moved all project files to new folder
  ([`4ad116f`](https://github.com/Baharis/hikari/commit/4ad116f5fc71a8dc57c91fcecbcd94d1fd19b0ea))

- Moved all relevant scripts to taskmasters
  ([`4b579c5`](https://github.com/Baharis/hikari/commit/4b579c5039e3783c9f3ad3a7ba45d322902d0bb8))

- Moved completeness maps to taskmasters, removed some old scripts
  ([`c6faace`](https://github.com/Baharis/hikari/commit/c6faacecc423e7809c2d9071ee5320c90c0885fe))

- Moved HklFrame.Crystal to BaseFrame Made HklFrame a sub-class of BaseFrame
  ([`6845870`](https://github.com/Baharis/hikari/commit/68458704c466525c3371d8cd43aa0aa1015cdebf))

- Moved kesshou/scripts/crysalis to separate project and removed from kesshou as it does not fit
  central theme. Fixed resolution having to be given in A instead of A-1 radius in completeness_map
  script.
  ([`c3226a3`](https://github.com/Baharis/hikari/commit/c3226a306cd3080b04d7a6a3cbd2f85e1af881a2))

- Moved kesshou/scripts/png to separate project and removed from kesshou as it does not fit central
  theme
  ([`3e323bd`](https://github.com/Baharis/hikari/commit/3e323bd0b53fa128f19b16f39b908b3b840fbe49))

- Moved scripts to common module Added first docstrings to the strings
  ([`0c2f785`](https://github.com/Baharis/hikari/commit/0c2f7855062c91a8f8cbcfd5606a5c6f3f12ae9f))

- Optimised imports, basic documentation all code in hkl dataframe. Renamed data to table.
  ([`ffbcfaf`](https://github.com/Baharis/hikari/commit/ffbcfaf0de88025acb94cefb8e71c7698dc712a2))

- Pdf file renamed
  ([`8310efe`](https://github.com/Baharis/hikari/commit/8310efe2458c004c8c4c9334e315e4649f23ec3b))

- Prepared tutorial in readme
  ([`81216f5`](https://github.com/Baharis/hikari/commit/81216f5f475fb3541b99c283bc1cd8dcbe3069bc))

- Provided initial documentation to HklFrame and BaseFrame.
  ([`e30dd22`](https://github.com/Baharis/hikari/commit/e30dd22185bea46f94bc2d7bfe970ee5fa761dc4))

- Read dat introduced Read hkl reintroduced
  ([`fb6ebe2`](https://github.com/Baharis/hikari/commit/fb6ebe29cfa31d20fdd8e747189561db9c327c3f))

- Read dat introduced Read hkl reintroduced
  ([`d77b735`](https://github.com/Baharis/hikari/commit/d77b735768a539ff7953e73f651e75e4a3bab023))

- Read dat introduced Read hkl reintroduced
  ([`711a6dd`](https://github.com/Baharis/hikari/commit/711a6dd214e3c2ce21017aeb4ab8914f8e4b325a))

- Recalculating I and si to F and sf ready
  ([`ff3ca02`](https://github.com/Baharis/hikari/commit/ff3ca02f9fd4b41c1445b2d4d980de737973153f))

- Refactored all instances of "kesshou" to "hikari"
  ([`692049c`](https://github.com/Baharis/hikari/commit/692049c47707b0b29424a45bd2cd6e7c401095da))

- Removed all pycache files
  ([`b5a92e2`](https://github.com/Baharis/hikari/commit/b5a92e2ff3e4b626d79eecb5e2b38f62fd0f2e30))

- Removed redundant matmul functionality of SymmOps Moved potential symm ops matrices to docs
  ([`ec5496c`](https://github.com/Baharis/hikari/commit/ec5496c08765df116900cc1f228923dd2b33c9ac))

- Removed some redundant or unused code from .hkl file
  ([`94277be`](https://github.com/Baharis/hikari/commit/94277bebe988c281829094b65c01b5d36edb5554))

- Removed test files
  ([`a466f89`](https://github.com/Baharis/hikari/commit/a466f899d74d3cd69d8904ba8122244796b34d8a))

- Reshkl moved to HklArtise, removed 2D draw
  ([`355e0d6`](https://github.com/Baharis/hikari/commit/355e0d6ad5fa8e7fb691c54548809b20439f0021))

- Resolution calculation introduced
  ([`3dbd99a`](https://github.com/Baharis/hikari/commit/3dbd99aed2d77fbed414a3f043b24c83299fe544))

- Resolution calculation moved to "place" of method of hkl dataframe
  ([`230efd9`](https://github.com/Baharis/hikari/commit/230efd931a5d94db2a603c86d24568c1569ec01d))

- Restucturised for dataframes
  ([`494a5b8`](https://github.com/Baharis/hikari/commit/494a5b81e38fc11b1dd7a9322634de9c39f1c848))

- Rewritten to *Frames
  ([`eaf2497`](https://github.com/Baharis/hikari/commit/eaf2497f9f649d23931d7cbe4bed887563c1e25c))

- Rewritten to *Frames
  ([`04a582d`](https://github.com/Baharis/hikari/commit/04a582d82ac69e06901aa01fea70a8c3592c8b96))

- Script for checking completeness some degrees around vector
  ([`f60901f`](https://github.com/Baharis/hikari/commit/f60901f8747be069bfe7dd496765a644b77b0293))

- Scripts added
  ([`abf12a4`](https://github.com/Baharis/hikari/commit/abf12a47a978d79cce503f0871b33cb5a7117f0b))

- Scripts moved to main directory
  ([`4484a00`](https://github.com/Baharis/hikari/commit/4484a0062ac1d37de76cee3828c63fc4a1808080))

- Sg generation now relies only on json definitions
  ([`1128c9d`](https://github.com/Baharis/hikari/commit/1128c9da951e5474affae1643ece0e7dc75ec7a9))

- Some additions I already forgot about
  ([`213048a`](https://github.com/Baharis/hikari/commit/213048aaa2677b81a25abf2c57406d1cef93d8d4))

- Some files renamed, increased readability, added manual-type files
  ([`2d79c5e`](https://github.com/Baharis/hikari/commit/2d79c5e671f0cec2829a315796fa63372a51d6cb))

- Sphinx added. Dataframes moved to scripts A few minor methods removed from HklFrame
  ([`69eb346`](https://github.com/Baharis/hikari/commit/69eb346f6156c85f90b999b0afaf797fe671a625))

- Split HklIo into HklReader and HklWriter
  ([`d8f67be`](https://github.com/Baharis/hikari/commit/d8f67be0461d82d588ec449bd359bc3fe57e8b95))

- Superclassed pointgroup using group, added parity functions
  ([`3c70203`](https://github.com/Baharis/hikari/commit/3c7020398cf44fa98191272d1258165626e06dd2))

- Tools for merging RGB-layers, alt #14&#31 SG definitions, modified dac-stats script
  ([`445fbc0`](https://github.com/Baharis/hikari/commit/445fbc0993efb5eb3802305124d25c5eb23ea251))

- Transformation and translation are now stored internally using ints to avoid approximations
  ([`fbf0f27`](https://github.com/Baharis/hikari/commit/fbf0f275c4db220f0aab6060f890a08ef229aea6))

- Updated requirements, removed drop_zero, removed find_reflection, optimised trim
  ([`be82419`](https://github.com/Baharis/hikari/commit/be82419ab8df96d7991b4b0a20d0235868d49b75))

- Vectorised 'calculate_uncertainty' method
  ([`a65fae5`](https://github.com/Baharis/hikari/commit/a65fae5cc45e408dc77f98ec52e4b419a92763f0))

- Vectorised 'dac' method
  ([`172f95a`](https://github.com/Baharis/hikari/commit/172f95a790497fb5c9b114669161847c1126d715))

- Vectorised 'generate_ball' method and renamed it to 'make_ball'
  ([`4be0bdf`](https://github.com/Baharis/hikari/commit/4be0bdf2d0029c79078024c3c6b81c364afb19b4))

- Vectorised 'place' method
  ([`6f5d4fd`](https://github.com/Baharis/hikari/commit/6f5d4fde0da60c1719482b11736ed1d19362aa1a))
