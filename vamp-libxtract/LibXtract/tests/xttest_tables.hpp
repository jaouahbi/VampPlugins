

// We need a pre-computed uniform distribution because tests need to be reproducible
// If we use "real" randomness then in certain cases tests will pass or fail depending on the RNG output for that specifc run
const double xttest_noise1024[1024] =
{
    -0.926215445821,-0.098017781796,0.109966880378,0.820549466690,-0.546435658457,-0.504663461036,0.793264509872,-0.634431900219,
    -0.152810838178,-0.503402913010,0.915663440593,-0.973438530483,0.158277372886,-0.041142932741,-0.327499969152,-0.256058057458,
    -0.634202010132,0.899549568118,0.908553902996,-0.732737332942,-0.502742283623,-0.004679979574,-0.455615446359,-0.699428826908,
    -0.471109507809,0.225955961633,0.463281528132,-0.121305692385,0.996208108905,-0.555233225774,-0.476810920977,0.980720733532,
    0.634634558284,-0.433651049607,0.564862557910,0.663526437942,-0.510615645620,-0.470264920381,0.542477393583,0.090342163732,
    -0.513409624565,-0.470838307845,-0.231490101753,-0.257177708854,-0.163867239218,-0.083145959415,-0.022323519875,0.229480646721,
    0.719916798685,0.024505892335,0.614007860190,0.352966930486,-0.663065348312,0.392548137668,0.776646838692,0.330756817911,
    -0.014892953556,0.378561502033,-0.463619260706,0.592079595980,0.070040783872,0.381476877464,-0.963381037295,0.895344541414,
    0.404249720331,-0.398873011267,0.452715168815,-0.813986729784,0.041195962800,0.803346069103,-0.953649412358,-0.109967652169,
    0.204966320073,-0.418625445237,0.339055604990,-0.980214827615,0.926327979049,0.638821764258,0.907810717787,0.081812727135,
    -0.320151906833,-0.778897475493,-0.698424031828,0.727345492782,-0.196055536723,0.123627253484,-0.901948700487,-0.859950094964,
    0.503337679604,-0.116340838742,-0.824027239600,0.284649443573,-0.624950288448,-0.192645040549,-0.675992602199,0.294178953210,
    0.967850433521,-0.625393243657,-0.813689226568,-0.576785475603,-0.488479760627,0.082945250462,-0.530351262068,0.684113017039,
    -0.383119845566,0.025359700098,-0.068407749257,-0.748119766699,-0.368902342547,-0.197390778219,0.791592633974,-0.753477687297,
    -0.976229986793,-0.454210604959,0.958680734329,0.563726744666,0.242113130085,-0.123640315607,0.681353911710,0.311449093127,
    -0.315951038977,-0.446555748418,0.427263443459,0.330388637867,0.630965200922,0.979410892225,-0.226522167009,0.614949799667,
    0.922541849498,-0.795865193626,-0.469458603316,-0.847649515413,0.514317987720,0.335073571573,-0.153357162068,0.558567738666,
    -0.845079920390,-0.271874938183,-0.117099332867,0.275992308117,0.635742752118,0.163833609389,0.815861330886,0.050066826416,
    0.540996092040,-0.653513525382,-0.838619646533,-0.178367993156,-0.250437501384,-0.247653169528,0.623550291388,0.757995457506,
    0.646119793178,-0.080913995679,0.842345790373,-0.729131055201,0.369358109245,-0.740655331485,0.059613604626,0.868112788883,
    -0.854860769232,-0.719223578146,-0.736130035541,-0.924101598760,0.610395578976,0.946877893985,0.862259027246,0.357827553602,
    -0.961871397142,-0.897688462729,0.080600617602,0.641828798083,0.753076988801,-0.702491476631,0.601960803921,-0.242838518948,
    0.336049684297,0.692679395857,0.159634505915,-0.456103932689,-0.957801058149,0.158817635215,0.351750704253,0.440394618637,
    -0.343282558350,0.930456909413,-0.846058064429,-0.055446990852,0.212760806057,-0.139278147316,0.912814658771,0.143949659097,
    -0.806101433055,0.156287960703,-0.297319585426,0.707692489839,-0.727910176913,-0.035567449444,-0.306770253998,-0.412909394567,
    0.392850583762,-0.559614011731,-0.637401732914,0.229656399984,0.529401155867,-0.470410697936,-0.561330551502,0.243352992522,
    -0.345006926090,0.964777300037,-0.805162074863,0.043667972855,0.439249964177,0.789413816974,-0.852997969065,0.518113161501,
    -0.165364373540,0.632756835592,-0.717100938299,0.535536208942,-0.463983079269,0.437931371294,0.562503125575,0.793919465490,
    -0.403799785247,0.676961766086,-0.913616212637,0.590621351423,0.499343253508,-0.060981568241,0.772854379487,0.115814213644,
    -0.107137853754,-0.714310148527,0.338253927463,-0.340002031512,0.520615298021,0.386603894911,0.942877998750,-0.292295809169,
    0.169905364136,0.441407421789,-0.779780504695,-0.006443793096,-0.678883416320,0.765883077821,0.518030341539,0.705966052440,
    -0.450072806261,0.964733410450,0.125038295245,-0.690182090950,-0.565168129118,0.697156632599,0.497034841775,0.559372302503,
    -0.022135690731,0.884374779629,0.901081352511,-0.847110236457,-0.369796128830,0.000199076543,0.078078885868,0.220945499443,
    -0.749324797417,-0.196705487064,0.107855039504,0.387673025332,-0.007095171465,0.733044715980,-0.809756254220,0.539619219816,
    0.240125748564,0.892045154042,0.017049752868,-0.309879415790,0.677172046601,-0.357756888620,0.991358101988,0.971924716553,
    -0.422606886296,-0.028367312955,0.929031599699,0.342199730409,-0.852903759921,-0.155187142362,-0.505436378144,-0.013896062967,
    -0.630723729288,-0.012607424742,0.876074261564,-0.728224414316,-0.886378210338,0.263522034583,-0.716836001474,0.547440883857,
    -0.690236737462,-0.520449921476,-0.894227379197,-0.745297883684,-0.958101360366,0.077264514862,-0.640555957921,-0.688986488087,
    -0.771073386873,-0.841479545101,0.986318027395,-0.110689706388,-0.195029194763,0.703707848121,0.239389062632,-0.751743344589,
    0.215011039675,0.693668755759,-0.051212517905,0.312129514890,-0.825162844299,0.503614112572,0.984839130175,-0.303865545206,
    -0.483821212075,0.302609971108,-0.071350900248,-0.080996566660,0.081382255829,0.921721281592,-0.303696941426,-0.908226906896,
    0.664314443467,-0.905973674417,-0.782782287758,0.103366601606,0.982117842415,-0.950506545707,0.021375723278,-0.238653323131,
    -0.288783883049,0.188287033174,0.163551580525,0.916603284248,0.592898122672,0.961625677318,0.093110577135,0.831715141771,
    -0.947691479006,-0.459438816071,0.186945725769,0.629659407162,-0.979292747824,0.418629567443,-0.287504368876,-0.225442707525,
    0.882794564794,0.492696327363,0.660526969994,-0.357532941914,0.122312285023,0.954270249908,0.471984162653,0.886974182514,
    -0.381433640044,0.127962992957,0.150559859842,0.286761746099,0.773062980574,0.009398146896,0.662805145398,0.802357002476,
    0.853051427056,-0.305848892703,0.744636870331,-0.044608459006,-0.292149419020,-0.672439028738,-0.559598630638,0.999127733464,
    0.543959368977,-0.107820084972,-0.710322348211,0.346335646111,-0.855913709820,-0.231566227380,0.696632518624,-0.391462816129,
    0.777916503711,0.523436403372,0.237756208218,0.854658504277,0.462588376310,-0.299590886473,-0.159631079694,-0.159286851275,
    -0.044380556306,0.831899263002,0.698436006608,-0.194771868090,0.526302641386,0.685551231118,0.190047000297,-0.391902389054,
    -0.071123030643,0.850144821166,-0.387388249837,0.179164096950,0.361511202007,-0.099387085185,-0.642351304895,0.943523699850,
    -0.952170038562,0.508566721904,0.090150414214,0.237559777496,-0.371678129187,0.457791216910,-0.676845530765,0.174989308392,
    0.779936916704,0.051320370414,0.118445866993,-0.918192140143,0.859285293614,0.040847095678,-0.316802266088,0.104997850078,
    0.486341281922,-0.684787856047,-0.819177025290,-0.421845079873,-0.179264543684,0.527338462798,-0.112654257152,-0.989722421149,
    0.704943983422,-0.705451942021,0.401688733519,0.627849956254,-0.603362149219,0.137513085839,0.809283744577,-0.456727896753,
    -0.695517805034,0.294281103914,-0.637668831508,0.697367439498,-0.360873542819,0.245165795706,0.583761200949,0.296360915835,
    0.974806392886,0.232077959736,0.562764815024,-0.350069366852,0.171650255017,-0.335783903344,-0.321101974152,-0.729900734558,
    0.080044248719,0.644806087371,0.755415289199,0.396520515713,0.215029369671,0.813583578774,0.805375483064,-0.628351026800,
    0.950890996049,-0.217263535262,0.636601559367,0.104576946835,0.660218109850,0.674945260898,-0.080514742927,-0.012984100428,
    -0.785841932353,-0.137510667518,0.561589931649,0.781099323878,-0.029917789277,-0.872116015920,0.200068837912,0.024966664994,
    0.902272952533,-0.052846681611,0.395558460686,-0.912735522798,-0.142049476590,0.448044177366,-0.844382585678,0.879224783326,
    0.530682080800,-0.765404470266,0.526193193800,0.975747760183,0.899995352230,0.902229779202,-0.954115343092,0.100205296330,
    0.321177325844,0.110160914120,0.553046007408,0.704007783151,0.180093550752,-0.287556327047,0.320020837205,0.768965659958,
    -0.727384378126,0.835627097539,-0.363401469452,0.017239462625,0.338445771225,0.804999063597,0.491333683543,-0.382225772449,
    -0.139502073414,-0.219212399014,0.351543996077,-0.563571744985,-0.219360956815,-0.205192899063,0.634182689920,-0.435653811180,
    0.511889081812,0.106907451790,-0.152891430287,0.535411742181,0.865943749704,-0.884904981386,-0.452288140821,-0.180277296455,
    0.376910707932,-0.272793359851,0.894713770731,0.137607596641,0.505992486726,-0.695215646021,0.092418893892,-0.151536037658,
    -0.192280464391,-0.392318391315,-0.958923699038,0.450712320444,-0.990308235160,0.216031491577,-0.230399036144,0.325208353824,
    -0.131376442598,0.064223134754,-0.114112095069,0.559469934545,-0.094985212687,0.241575743405,-0.765110789890,0.525377582940,
    -0.110442579905,0.048886241227,0.453061747919,0.628324948234,0.061482120874,0.224216934420,-0.939834826328,0.083644032354,
    -0.077652651210,0.555337137678,-0.394845550230,0.926236412126,-0.227370305771,-0.276010098023,-0.909639516648,-0.655294507648,
    0.011439017637,-0.578912544773,-0.644236101516,0.013695945253,-0.017039548370,-0.239735559464,-0.789253266911,0.595443252403,
    0.705434906499,-0.409675565509,-0.155773560956,0.249490749536,-0.865745129895,-0.508055655762,0.583321268027,-0.817442297968,
    0.785291418495,-0.631002474398,0.491668229790,-0.075635239985,0.154131705807,-0.395183800493,-0.502275086479,-0.434708199075,
    0.123508432439,0.424565637114,0.564542812236,0.025784148531,0.976096746409,0.302736166686,-0.222625565782,0.205221138361,
    0.457932187633,0.169571765235,0.074500980661,-0.473512416194,-0.393208474553,-0.254238110737,-0.864602911402,-0.142201339153,
    -0.756447974523,0.958511251703,-0.553415877672,0.538713562583,-0.867953201277,0.258449408759,0.741997815725,-0.661583718467,
    -0.723028891040,0.116828769370,-0.450503654431,-0.966097125540,0.008718913612,-0.092277402764,-0.764685123178,-0.731776097221,
    0.188613372355,0.883880505874,0.501883506386,-0.510430726574,0.145045388882,-0.314989741915,0.198638248758,0.811702726099,
    -0.576229396322,-0.992168308358,0.010549605415,0.595176315765,-0.879916828582,-0.637384034311,-0.702623787644,-0.204115427445,
    -0.517120043959,-0.538593051574,0.481679443514,-0.543318746251,-0.597206481200,-0.351050092458,-0.593859575937,0.056037203711,
    -0.573506161501,0.192140154822,0.011416899192,0.177927272460,0.465092313177,0.962796252347,-0.327035417453,-0.689752729870,
    0.112788254775,0.553662726777,-0.347898153871,0.456436939766,-0.444350262761,0.842359584978,0.714299584899,-0.504617744416,
    -0.831502484889,0.570080898656,-0.073833523013,-0.604405918303,0.424747003921,0.310899611143,-0.764538095956,-0.416069816420,
    -0.355637687429,0.683954161272,0.911679770669,0.810664288224,0.990704889956,-0.641820612795,-0.664507255969,0.595613010707,
    -0.441659543489,-0.764871131739,-0.480955508046,-0.245017682899,-0.472803285492,0.129488866380,-0.456430610042,0.177972667187,
    -0.825477332513,0.025501833836,-0.133990018280,-0.770695398197,0.360192345825,0.792833757494,0.727201589852,0.521851799015,
    -0.671175573048,0.094261308772,0.279623189102,-0.030505536081,-0.017473004872,0.814894741434,-0.037787074232,-0.737475680834,
    0.664669262135,0.935077024437,0.104794700721,-0.097255238986,0.876209570749,0.950096908360,-0.868175242112,0.890058279453,
    0.324230527533,-0.605292658820,0.813575711806,0.640223535245,0.637834512602,-0.602396416293,0.702816913560,-0.406529356853,
    0.623918859553,0.471089933204,0.515379793395,0.959717164474,-0.894039690467,-0.609542778443,-0.451674298593,0.593383463790,
    0.570161280470,0.834312325312,-0.419466477336,0.647746453803,-0.577620190509,-0.054118193455,0.943808619129,0.557665303748,
    -0.106358104685,-0.597506351807,0.115787379373,-0.484386400581,-0.818955346506,-0.494904198883,0.900480444997,0.478306401639,
    0.532769709262,0.252769924896,-0.116577199708,-0.663499593221,-0.074125571939,0.229396860221,0.454739683436,0.715923268404,
    0.182655699039,-0.614598816017,-0.120158637688,0.227931764220,-0.541801481088,0.611728012935,-0.126609765720,-0.868739359465,
    -0.966686219368,0.244469143429,0.687962675120,0.578015366460,0.756819082474,0.207157081837,-0.611160741663,0.810133655907,
    -0.635041257469,-0.658538203971,0.487217191361,0.820352316884,-0.734683954105,0.600506273787,0.766671907232,-0.459316692596,
    0.357927947582,0.362353920762,0.929083260966,0.587683181111,-0.044784075488,0.709828798177,-0.849638396016,0.146387642226,
    0.101904507907,-0.581958485709,0.047882316958,-0.674891836456,-0.657701621662,0.449383087118,0.794652820542,-0.896764726852,
    -0.067837810813,-0.981097963726,0.635982815929,0.941128988628,0.175298967024,0.209466964232,0.797476082755,-0.256787620113,
    0.611580850331,0.800908538062,0.175184771023,-0.179378464982,-0.310512828643,0.043141609223,0.681149363550,-0.793087340036,
    -0.529126905197,-0.303324113583,0.055699083465,-0.380875661435,-0.645176444724,0.105365790845,-0.850734299121,0.384352816320,
    -0.448460750911,0.760241283573,-0.229023537369,0.431783936089,0.651391885657,-0.842881821715,-0.261095603402,0.802422348439,
    -0.483778819408,0.546412378552,-0.697745515579,-0.816307621646,0.860617524329,0.029774928752,-0.590532158614,-0.744522469238,
    0.211133380072,-0.454349305040,-0.982690751931,-0.835589513313,0.102974052993,0.794721275333,-0.807174031783,0.147439753030,
    -0.782129412275,0.131675788312,0.139092582048,-0.608747172309,-0.396084822597,0.159022790772,0.920212202166,-0.122311433720,
    0.160767085579,0.023735403496,0.560427865499,0.535353471070,0.765260438242,-0.446169272200,-0.538681366997,-0.197100455708,
    0.447626039654,-0.253307189810,-0.478292044515,-0.850853633835,0.859721788549,-0.118688796712,-0.361553269969,-0.021676787636,
    0.901751985855,0.306400058222,0.128373480250,-0.046004164430,0.526728627132,-0.608455444127,-0.881979243730,-0.372620944720,
    -0.657508461124,0.732804249959,0.581334739811,0.777764571675,0.843461999903,0.584723829362,-0.365025662958,-0.760186457886,
    0.837157659372,-0.060462860641,-0.422567815531,0.117466888651,0.882441887322,0.737637086738,-0.462941365372,0.150449364282,
    0.266552755794,-0.063535223507,-0.628314426898,-0.341507243988,0.476528395415,-0.127086605866,0.419293439132,-0.432213558536,
    -0.106550569079,0.210175348183,-0.678715074359,-0.899176692011,-0.130227874266,0.144548908730,0.841204926662,0.889456666706,
    0.581527425053,0.463261416608,0.447525061523,0.278247400874,0.983593223205,0.642710053683,0.238582318740,-0.282250205731,
    -0.791158602728,-0.462884996568,0.853019622525,-0.635694418685,-0.856674102673,-0.775044025526,-0.835518826278,-0.692627215976,
    0.985006310183,0.963524706440,0.529960969668,0.113165633509,0.758963583933,-0.453414896524,-0.138587961787,-0.711175920330,
    0.847713775861,-0.703713632218,-0.517598751164,0.947025826371,-0.485635710655,-0.351566647084,0.443855276858,0.493667233018,
    -0.523834569548,0.869953567440,-0.359673750456,0.219182349556,-0.297940802507,-0.497006232259,0.043814088923,-0.070101336133,
    0.286655727890,0.526934875827,-0.786852733839,-0.964649559846,0.570214953890,0.517424334124,-0.909379614720,-0.946028492056,
    0.966222106608,0.782350997081,0.771620835482,-0.826545737777,-0.554431775276,0.050357634271,0.198384108950,-0.031026258856,
    0.673880678005,0.640256311558,0.522082402766,0.394933672318,0.496403687157,0.694107283834,-0.349141222251,-0.506073162822,
    0.046303361784,-0.304125595493,0.175011905085,-0.279951886878,-0.842643261884,-0.543541406683,0.153043694999,0.603749866585,
    0.776056479698,-0.107985838547,-0.002847460038,0.659511743503,0.055756607387,-0.818144588874,-0.231981372610,0.229869051308,
    0.683206080637,0.303546289083,-0.653814889780,0.535763370295,-0.020375278904,0.792879778451,0.585093338499,-0.308757313343,
    -0.257467533251,0.701683251995,-0.391768340938,0.522651455997,-0.288019201002,-0.984653746978,-0.739779482875,-0.155129086557,
    0.674317270393,-0.518343751954,-0.733193704902,0.764046753044,0.358910044856,0.872947514948,-0.554140649693,-0.234887097514,
    -0.749812536925,-0.563867423759,0.703964695518,0.261991771960,-0.525598087382,0.435382571368,0.221018111661,0.870558675422,
    -0.198951908611,-0.777248756421,0.835998664352,-0.724856869417,-0.747737442734,-0.968655476872,0.968473730666,0.230038567716
};