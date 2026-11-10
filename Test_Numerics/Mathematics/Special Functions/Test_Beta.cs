using System;
using Microsoft.VisualStudio.TestTools.UnitTesting;
using Numerics.Mathematics.SpecialFunctions;

namespace Mathematics.SpecialFunctions
{
    /// <summary>
    /// Unit tests for the Beta class
    /// </summary>
    /// <remarks>
    ///      <b> Authors: </b>
    /// <list type="bullet">
    /// <item><description>
    ///     Haden Smith, USACE Risk Management Center, cole.h.smith@usace.army.mil
    /// </description></item>
    /// <item><description>
    ///     Sadie Niblett, USACE Risk Management Center, sadie.s.niblett@usace.army.mil
    /// </description></item>
    /// </list>
    /// <b> Description: </b>
    /// Some of these special functions are tested directly with MS Excel and with values provided by:
    /// "Numerical Recipes, Routines and Examples in Basic", J.C. Sprott, Cambridge University Press, 1991.
    /// Others must be tested indirectly through the continuous distribution classes.
    /// </remarks>
    [TestClass]
    public class Test_Beta
    {
        /// <summary>
        /// Test the Beta Function
        /// </summary>
        [TestMethod]
        public void Test_Function()
        {
            var testX = new double[] { 1d, 0.2d, 1d, 0.4d, 1d, 0.6d, 0.8d, 6d, 6d, 6d, 6d, 6d, 7d, 5d, 4d, 3d, 2d };
            var testY = new double[] { 1d, 1d, 0.2d, 1d, 0.4d, 1d, 1d, 6d, 5d, 4d, 3d, 2d, 7d, 5d, 4d, 3d, 2d };
            var testValid = new double[] { 1d, 5d, 5d, 2.5d, 2.5d, 1.666667d, 1.25d, 0.0003607504d, 0.0007936508d, 0.001984127d, 0.005952381d, 0.0238095d, 0.00008325008d, 0.001587302d, 0.007142857d, 0.03333333d, 0.1666667d };
            var testResults = new double[testValid.Length];

            for (int i = 0; i < testValid.Length; i++)
            {
                testResults[i] = Beta.Function(testX[i], testY[i]);
                Assert.AreEqual(testResults[i], testValid[i], 1E-4);
            }
        }

        /// <summary>
        /// The log beta function matches 400-digit reference values on every branch of algorithm 708, including large,
        /// nearly equal shapes and strongly asymmetric shapes, where ln Γ(a) + ln Γ(b) - ln Γ(a + b) loses most of its
        /// digits.
        /// </summary>
        /// <remarks>
        /// Reference values computed with Python mpmath (1.3.0, and 1.4.1 for the a &lt; 1, 2 &lt; b &lt; 8 points) at 400 significant digits as
        /// loggamma(a) + loggamma(b) - loggamma(a + b). The tolerance is 1E-15 relative to the larger of 1 and the result,
        /// about five ulps; at (0.5, 1E6) the difference of ln Γ values is in error by about 1E-9, and at (0.5, 1E15)
        /// it returns -16 instead of -16.697. The points with a &lt; 1 and 2 &lt; b &lt; 8 cover the branch that reduces b
        /// to (1, 2], where ln Γ(b) - ln Γ(a + b) taken directly loses up to about 25 ulps.
        /// </remarks>
        [TestMethod]
        public void Test_LogFunction_ReferenceValues()
        {
            // { a, b, ln B(a, b) }
            var table = new double[,]
            {
                { 0.5, 0.5, 1.1447298858494002 },
                { 2, 3, -2.4849066497880004 },
                { 0.001, 0.002, 7.313217104417236 },
                { 1E-10, 3.5, 23.02585092977242 },
                { 0.5, 1E6, -6.335390211057437 },
                { 0.5, 1E15, -16.69702325453064 },
                { 3.7, 1E7, -58.208882081380345 },
                { 7.9, 7.9, -10.703832437848734 },
                { 1.5, 1500.25, -11.091112700256792 },
                { 5.5, 2500, -39.07938579535403 },
                { 1000, 1000, -1388.4826016359023 },
                { 1E5, 1E5, -138633.92706134808 },
                { 8, 1E10, -175.68164608125824 },
                { 1E8, 0.25, -3.317147660352514 },
                { 2.5E7, 3E7, -37895515.41140297 },
                { 12.5, 1E300, -8615.959751215734 },
                { 1E-300, 1E-300, 691.4686750787737 },
                { 0.9219, 7.5939, -1.8140586851054754 },
                { 0.7052453993501302, 7.681687846760826, -1.16997303546981 },
                { 0.8739311446570729, 7.6339892317246765, -1.6825461302652165 },
                { 0.5, 2.5, 0.16390063283767395 },
                { 0.25, 3.75, 0.9830786340634395 },
                { 0.1, 5.5, 2.0906131618193924 },
                { 0.999, 2.001, -0.6921460352922537 },
                { 0.01, 7.99, 4.5793301230650885 },
                { 0.75, 4.25, -0.8603159514662788 },
            };
            for (int i = 0; i < table.GetLength(0); i++)
            {
                double a = table[i, 0], b = table[i, 1], expected = table[i, 2];
                double tolerance = 1E-15 * Math.Max(1d, Math.Abs(expected));
                Assert.AreEqual(expected, Beta.LogFunction(a, b), tolerance, $"a = {a:R}, b = {b:R}");
                Assert.AreEqual(Beta.LogFunction(a, b), Beta.LogFunction(b, a), 0d, $"symmetry at a = {a:R}, b = {b:R}");
            }
        }

        /// <summary>
        /// ln B(a, b) on 20 random points in each region of the documentation table, with a the smaller and b the larger
        /// shape from 1E-10 to 1E300, matches mpmath (ln Γ(a) + ln Γ(b) - ln Γ(a + b) at 400 digits) to the table bound,
        /// relative to max(1, |ln B|). The largest error measured on these points is 5.6E-16, for a &lt; 1 and 2 &lt; b &lt; 8.
        /// </summary>
        [TestMethod]
        public void Test_LogFunction_RegionBounds_MpmathOracle()
        {
            var cases = new (double a, double b, double expected, double bound)[]
            {
                // a < 1, b ≤ 2: at most 2.8E-16
                (4.425629092419659e-05, 0.7606045372139609, 10.025534358838303, 2.8E-16),
                (0.011946652803248529, 0.29977015759671194, 4.4615478856309236, 2.8E-16),
                (0.00012027142827833866, 0.48769125330360835, 9.0259336358644298, 2.8E-16),
                (1.4978780409531917e-08, 0.14272752393865532, 18.016631378510976, 2.8E-16),
                (3.727424940073971e-05, 0.12926033372749338, 10.197488904299254, 2.8E-16),
                (1.3555011306315875e-08, 0.15004012419346058, 18.116509607210041, 2.8E-16),
                (0.37385546172693257, 0.7167646150370708, 1.1528730468137041, 2.8E-16),
                (4.719105075154764e-09, 0.4717616429147073, 19.171646665280434, 2.8E-16),
                (7.648832099615353e-08, 1.564553621052137, 16.386127723056067, 2.8E-16),
                (1.825127373381545e-06, 1.2586099454411652, 13.213860122646039, 2.8E-16),
                (0.3131213025168426, 1.6231293768113675, 0.96704272531003488, 2.8E-16),
                (0.004772497573377027, 0.15105305366696206, 5.3749244276256789, 2.8E-16),
                (0.00033891006351785004, 1.2528762298936513, 7.9896561078718209, 2.8E-16),
                (2.094924543558192e-08, 1.6758379931076746, 17.681163192565666, 2.8E-16),
                (4.186254407536011e-10, 0.8528145548708301, 21.594044532068473, 2.8E-16),
                (1.550638205991614e-05, 0.8939886170670407, 11.074261808717191, 2.8E-16),
                (3.371777159048381e-09, 1.9986570886260466, 19.507825881313167, 2.8E-16),
                (2.279331976066434e-09, 1.0324474882477674, 19.899383429745349, 2.8E-16),
                (0.5420407382156988, 1.304753637518221, 0.44196546973082163, 2.8E-16),
                (0.0019227097913468954, 1.9130411388883548, 6.2522096828934881, 2.8E-16),
                // a < 1, 2 < b < 8: at most 6.7E-16
                (1.9165252651151692e-05, 7.966335862968, 10.862362069878742, 6.7E-16),
                (0.015087543336333189, 5.676035731677544, 4.1605131816627638, 6.7E-16),
                (1.1482146199255915e-08, 2.393553116359275, 18.282472498287298, 6.7E-16),
                (6.912085323550568e-06, 2.6384531393908093, 11.882229877983868, 6.7E-16),
                (4.064764464668389e-08, 3.9026599666099693, 17.018324871520313, 6.7E-16),
                (3.933151907324852e-07, 5.432595829148478, 14.748653679100510, 6.7E-16),
                (0.23255782929185093, 2.1355880002990713, 1.2312621482723967, 6.7E-16),
                (5.417012417340264e-07, 2.0486391225177822, 14.428550643554390, 6.7E-16),
                (1.0949505220642757e-10, 3.1488453327747896, 22.935141752851809, 6.7E-16),
                (0.012434520214830254, 5.906661700530784, 4.3592110302177618, 6.7E-16),
                (4.0908014212187275e-09, 6.851217434196128, 19.314524929650535, 6.7E-16),
                (0.0009433073014158004, 2.775861650917615, 6.9647914083262714, 6.7E-16),
                (0.725061489136191, 3.3709963929394497, -0.62110120420882331, 6.7E-16),
                (0.004057777365110821, 2.0157253181126604, 5.5030294355542786, 6.7E-16),
                (0.0003105783021295323, 5.542046730678701, 8.0763923438050611, 6.7E-16),
                (0.00018999409788543554, 4.014147911156684, 8.5681684907148483, 6.7E-16),
                (0.014032567738282188, 6.837889140296927, 4.2324936605599364, 6.7E-16),
                (0.13261539981241136, 2.3813848228142085, 1.8675803032536802, 6.7E-16),
                (2.493136750961895e-06, 3.2230925665149943, 12.901964950633447, 6.7E-16),
                (0.2907565241151271, 6.473686554277267, 0.60173281531492739, 6.7E-16),
                // 1 ≤ a ≤ 2, b ≤ 2: at most 4.5E-16
                (1.7190043777838073, 1.7578990922428612, -1.3497071871252764, 4.5E-16),
                (1.5667176390058104, 1.881774451071962, -1.3062613180397242, 4.5E-16),
                (1.3147050958857978, 1.6687275809900095, -0.89039270652868174, 4.5E-16),
                (1.6743686317019204, 1.7629883941863538, -1.3145346471105250, 4.5E-16),
                (1.7054645055434485, 1.7070627416446285, -1.2947248504243125, 4.5E-16),
                (1.3074368724926186, 1.931005976789951, -1.0610862462801236, 4.5E-16),
                (1.3077946626842416, 1.5736708344028418, -0.81165369548469779, 4.5E-16),
                (1.6765512147485673, 1.8371129557181503, -1.3765604836654943, 4.5E-16),
                (1.0348332928047514, 1.6138211758397498, -0.52439618095776086, 4.5E-16),
                (1.25901845663233, 1.6782389916254563, -0.83646474912647040, 4.5E-16),
                (1.694390100299742, 1.8480720640528774, -1.4016193462453086, 4.5E-16),
                (1.0127911607800102, 1.7840234044765992, -0.59681963215568623, 4.5E-16),
                (1.0197189389922512, 1.517824896634854, -0.44247497748643155, 4.5E-16),
                (1.6114411872893788, 1.7803492213681795, -1.2711937054417124, 4.5E-16),
                (1.3791101942069255, 1.9363835088642687, -1.1468528522434487, 4.5E-16),
                (1.4480688145298013, 1.8849117724279996, -1.1870959520902367, 4.5E-16),
                (1.6473599545495998, 1.9176872539035679, -1.4116736858810597, 4.5E-16),
                (1.0626327417618122, 1.4178545301747663, -0.42462306517293807, 4.5E-16),
                (1.0837509556889944, 1.6471004232744455, -0.60812526447496251, 4.5E-16),
                (1.404887737606367, 1.774021796521616, -1.0625584603367346, 4.5E-16),
                // 1 ≤ a < 8, 2 < b < 8: at most 4.1E-16
                (3.222213661730355, 7.610457530071265, -6.0542894278824253, 4.1E-16),
                (3.472660186626311, 4.040538248640432, -4.5462208559368653, 4.1E-16),
                (3.568781126519939, 7.9487819325233335, -6.6344101790282970, 4.1E-16),
                (2.4764376901854277, 6.675667010497509, -4.6832458290842141, 4.1E-16),
                (2.1405385369973473, 4.972608893714716, -3.5896391360258462, 4.1E-16),
                (4.311330125965102, 7.308974531751, -7.2198690172859550, 4.1E-16),
                (3.079888547024888, 3.458149057019411, -3.8076696027965555, 4.1E-16),
                (2.6606275368630374, 4.0040503151839175, -3.7594019462272934, 4.1E-16),
                (3.7300032167353985, 5.076678106812512, -5.4355726702845732, 4.1E-16),
                (2.5337370840329796, 3.750948756064843, -3.4838035457089496, 4.1E-16),
                (3.0097392575855655, 5.850778247693356, -5.0700781244448283, 4.1E-16),
                (2.3978373604140586, 7.030651461763216, -4.6802798183025667, 4.1E-16),
                (3.3597643760706264, 4.11755300072663, -4.4993127583124233, 4.1E-16),
                (3.3518311570876924, 7.678584082076615, -6.2504997786143885, 4.1E-16),
                (3.45137614959896, 4.369605281353499, -4.7439753376400455, 4.1E-16),
                (2.8943482759846737, 4.3858156931918675, -4.2152047099568062, 4.1E-16),
                (1.239123009320278, 5.812353609543418, -2.3010114834233839, 4.1E-16),
                (1.3084162693398946, 3.8901983874730472, -1.9354879793971596, 4.1E-16),
                (2.98951065793983, 3.111263004750596, -3.4786202529430286, 4.1E-16),
                (2.455783740817233, 4.851170241897103, -3.9507807230848035, 4.1E-16),
                // a < 8, b ≥ 8: at most 1.2E-15
                (5.650638322498572e-08, 9.734428599875051e+217, 16.688883832501940, 1.2E-15),
                (2.7808838846168915e-06, 9.530285170132884e+98, 12.792106346358622, 1.2E-15),
                (0.0011924632485075149, 1.711457750339112e+199, 6.1840023782714253, 1.2E-15),
                (6.68141675145515e-10, 3.30467307535333e+139, 21.126520661309016, 1.2E-15),
                (0.00015547926143552054, 1.0009926331996369e+208, 8.6944434430039070, 1.2E-15),
                (4.225188653371009e-09, 6.510834007858554e+72, 19.282201214616425, 1.2E-15),
                (5.407606640294096e-06, 8316052038610845.0, 12.127502610727785, 1.2E-15),
                (1.1381380324272365e-05, 158241484.21049455, 11.383310397063028, 1.2E-15),
                (0.0766607864641367, 9.057731292465986e+33, -3.4652476180814200, 1.2E-15),
                (1.1589912846278578e-09, 1.1876970486178249e+200, 20.575715257747881, 1.2E-15),
                (8.732554701519961e-05, 7.39350144577857e+87, 9.3281489184009363, 1.2E-15),
                (8.869228403635027e-10, 4.0833674963541626e+200, 20.843262716643558, 1.2E-15),
                (6.0291502496908475e-09, 4.963141422218253e+247, 18.926656314437786, 1.2E-15),
                (1.329788863089703e-08, 9.178919936946513e+260, 18.135652565646161, 1.2E-15),
                (2.2528883176852884e-10, 3.1849508298883696e+109, 22.213637783274923, 1.2E-15),
                (1.0191875955186441e-09, 1.8153124647730776e+165, 20.704259613574266, 1.2E-15),
                (0.02256525144986473, 7.037570438256299e+166, -4.8903937239135606, 1.2E-15),
                (4.680638114889795e-10, 1.3767209590656015e+272, 21.482416186400549, 1.2E-15),
                (9.577066047825932e-06, 3.581934785934263e+210, 11.551490601580053, 1.2E-15),
                (0.013453034825016674, 1.2767632099457602e+146, -0.22496039657007055, 1.2E-15),
                // a ≥ 8: at most 2.6E-16
                (5.1513173283645065e+191, 2.9278503648061483e+258, -7.9695112034277259E+193, 2.6E-16),
                (5.550585747626707e+176, 5.8918566927198506e+193, -2.2315360850642222E+178, 2.6E-16),
                (418818749831.7738, 1.7366870413034903e+153, -136990082334626.81, 2.6E-16),
                (7.915259786323023e+31, 3.0765666688773996e+66, -6.3833008307655250E+33, 2.6E-16),
                (5.2673902375280465e+140, 1.2374294230234936e+159, -2.2808125211084783E+142, 2.6E-16),
                (9.453084608893278e+205, 6.229260906156618e+220, -3.3200833538667521E+207, 2.6E-16),
                (6.06354001967649e+62, 1.067886095400101e+229, -2.3271569349958381E+65, 2.6E-16),
                (5278254182738227.0, 4.954178592860689e+26, -1.3863372618486854E+17, 2.6E-16),
                (4.724709323246599e+149, 2.2248495708974415e+208, -6.4303012179174986E+151, 2.6E-16),
                (8.324351134895498e+179, 2.1760439563023236e+196, -3.2300372047018671E+181, 2.6E-16),
                (1.6457329510180383e+104, 2.092350784928545e+124, -7.7829680523435774E+105, 2.6E-16),
                (3.3967829369309535e+20, 7.708435913818512e+239, -1.7190630050260353E+23, 2.6E-16),
                (1.4682305772147003e+31, 7.754692299731448e+158, -4.3326390052784715E+33, 2.6E-16),
                (4.86118301088611e+53, 1.4166958502143365e+271, -2.4390042480566543E+56, 2.6E-16),
                (2.028738378541322e+78, 3.620096866018765e+291, -9.9819955906276684E+80, 2.6E-16),
                (1.349420966084688e+192, 3.218312936392213e+257, -2.0448749202370398E+194, 2.6E-16),
                (5.372092145255996e+21, 3.272290475365136e+187, -2.0560790821786083E+24, 2.6E-16),
                (2.35148810963231e+78, 1.2968837489339554e+282, -1.1055104543240769E+81, 2.6E-16),
                (4.290099986177925e+48, 8.102995165288503e+84, -3.6263780628268084E+50, 2.6E-16),
                (3.9700859605934975e+220, 5.399486461570942e+269, -4.5312252810872636E+222, 2.6E-16),
            };
            foreach (var (a, b, expected, bound) in cases)
            {
                Assert.AreEqual(expected, Beta.LogFunction(a, b), bound * Math.Max(1d, Math.Abs(expected)), $"a = {a:R}, b = {b:R}");
            }
        }

        /// <summary>
        /// The log beta function returns NaN for NaN, negative infinity for a positive infinite argument, and
        /// ln|Γ(a)| + ln|Γ(b)| - ln|Γ(a + b)| for non-positive arguments; the beta function is its exponential.
        /// </summary>
        /// <remarks>
        /// scipy.special.betaln(-0.5, 2.0) = 1.3862943611198906 (ln 4) and betaln(inf, 2.0) = -inf.
        /// </remarks>
        [TestMethod]
        public void Test_LogFunction_SpecialValues()
        {
            Assert.IsTrue(double.IsNaN(Beta.LogFunction(double.NaN, 2d)));
            Assert.IsTrue(double.IsNaN(Beta.LogFunction(2d, double.NaN)));
            Assert.AreEqual(double.NegativeInfinity, Beta.LogFunction(double.PositiveInfinity, 2d));
            Assert.AreEqual(double.NegativeInfinity, Beta.LogFunction(0.5d, double.PositiveInfinity));
            Assert.AreEqual(1.3862943611198906, Beta.LogFunction(-0.5d, 2d), 1E-15);
            Assert.AreEqual(Gamma.LogGamma(-0.5d) + Gamma.LogGamma(2d) - Gamma.LogGamma(1.5d), Beta.LogFunction(-0.5d, 2d), 0d);
            foreach (var (a, b) in new[] { (0.3d, 0.7d), (2d, 3d), (7.5d, 12.25d), (40d, 60d) })
            {
                Assert.AreEqual(Math.Exp(Beta.LogFunction(a, b)), Beta.Function(a, b), 0d);
            }
        }

        /// <summary>
        /// Test the Incomplete Beta function Ix(a,b)
        /// </summary>
        /// <remarks>
        /// This code was copied from the Accord Math Library
        /// <para></para>
        /// References:
        /// <list type="bullet">
        /// <item><description>
        /// Accord Math Library, <see href = "http://accord-framework.net" />
        /// </description></item>
        /// </list>
        /// </remarks>
        [TestMethod]
        public void Test_Incomplete()
        {
            double a = 5;
            double b = 4;
            double x = 0.5;

            double actual = Beta.Incomplete(a, b, x);
            double expected = 0.36328125;

            Assert.AreEqual(expected, actual, 1e-6);
        }

        /// <summary>
        /// Test the continued fraction expansion #1 for the incomplete beta integral
        /// </summary>
        /// <remarks>
        /// This code was copied from the Accord Math Library.
        /// <para></para>
        /// References:
        /// <list type="bullet">
        /// <item><description>
        /// Accord Math Library, <see href = "http://accord-framework.net" />
        /// </description></item>
        /// </list>
        /// </remarks>
        [TestMethod()]
        public void Test_Incbf()
        {
            double actual = Beta.Incbcf(4, 2, 4.2);
            Assert.AreEqual(-0.23046874999999992, actual);
        }

        /// <summary>
        /// Test the continue fraction expansion #2 for incomplete beta integral
        /// </summary>
        /// <remarks>
        /// This code was copied from the Accord Math Library.
        /// <para></para>
        /// References:
        /// <list type="bullet">
        /// <item><description>
        /// Accord Math Library, <see href = "http://accord-framework.net" />
        /// </description></item>
        /// </list>
        /// </remarks>
        [TestMethod()]
        public void Test_Incbd()
        {
            double actual = Beta.Incbd(4, 2, 4.2);
            Assert.AreEqual(0.7375, actual);
        }

        /// <summary>
        /// Test the inverse of the incomplete beta integral
        /// </summary>
        /// <remarks>
        /// This code was copied from the Accord Math Library.
        /// <para></para>
        /// References:
        /// <list type="bullet">
        /// <item><description>
        /// Accord Math Library, <see href = "http://accord-framework.net" />
        /// </description></item>
        /// </list>
        /// </remarks>
        [TestMethod()]
        public void Test_IncompleteInverse()
        {
            double actual = Beta.IncompleteInverse(0.5, 0.6, 0.1);
            Assert.AreEqual(0.019145979066925722, actual);

        }

        /// <summary>
        /// Test the power series for the incomplete beta integral
        /// </summary>
        /// <remarks>
        /// This code was copied from the Accord Math Library.
        /// <para></para>
        /// References:
        /// <list type="bullet">
        /// <item><description>
        /// Accord Math Library, <see href = "http://accord-framework.net" />
        /// </description></item>
        /// </list>
        /// </remarks>
        [TestMethod()]
        public void Test_PowerSeries()
        {
            double actual = Beta.PowerSeries(4, 2, 4.2);
            Assert.AreEqual(-3671.801280000001, actual);
        }

        /// <summary>
        /// Test the incomplete beta function ratio
        /// </summary>
        /// <remarks>
        /// This function was test against values calculated using the "Ibeta( )" and "Cbeta( )" functions the R "zipfR" package
        /// <para></para>
        /// References:
        /// <para></para>
        /// Evert, Stefan and Baroni, Marco (2007). "zipfR: Word frequency distributions in
        /// R." In Proceedings of the 45th Annual Meeting of the Association for Computational
        /// Linguistics, Posters and Demonstrations Sessions, pages 29-32, Prague, Czech
        /// Republic. (R package version 0.6-70 of 2020-10-10)
        /// </remarks>
        [TestMethod()]
        public void Test_IncompleteRatio()
        {
            var testX = new double[] { 0.2d, 0.2d, 0.2d, 0.5d, 0.5d, 0.5d, 0.8d, 0.8d, 0.8d };
            var testP = new double[] { 2d, 3d, 4d, 2d, 3d, 4d, 2d, 3d, 4d};
            var testQ = new double[] { 4d, 6d, 9d, 4d, 6d, 9d, 4d, 6d, 9d};
            var testBeta = new double[testX.Length];
            for (int i = 0; i < testBeta.Length; i++)
            {
                testBeta[i] = Math.Log(Gamma.Function(testP[i]) * Gamma.Function(testQ[i]) / Gamma.Function(testP[i] + testQ[i]) ); 
            }
            var testValid = new double[] { 0.26272d, 0.2030822d, 0.2054311d, 0.8125d, 0.8554688d, 0.927002d, 0.99328d, 0.9987686d, 0.9999378d};
            var testResults = new double[testValid.Length];

            for(int i = 0; i < testResults.Length; i++)
            {
                testResults[i] = Beta.IncompleteRatio(testX[i], testP[i], testQ[i], testBeta[i]);
                Assert.AreEqual(testValid[i], testResults[i], 1e-4);
            }
        }

        /// <summary>
        /// One large and one small shape with x just below a / (a + b), where the continued fractions converge slowly:
        /// the asymptotic expansion matches mpmath (the continued fraction at 40 digits) to 1E-13 relative. The continued
        /// fractions alone were off by up to 1.5E-7 at these points.
        /// </summary>
        [TestMethod]
        public void Test_Incomplete_LargeShape_MpmathOracle()
        {
            var cases = new (double a, double b, double x, double expected)[]
            {
                (3217.498503960504, 0.052249964750908894, 0.9999836000654025, 0.12044471354300114),
                (1167.80589183603, 0.025166336558533956, 0.9999606071761743, 0.062560438832504961),
                (2305.3599894211384, 0.01839213353908136, 0.9999773323227816, 0.043888529386076856),
                (603.8638489638837, 0.021083286313687603, 0.9998821141221607, 0.044349845868886753),
            };
            foreach (var (a, b, x, expected) in cases)
            {
                Assert.AreEqual(expected, Beta.Incomplete(a, b, x), expected * 1E-13, $"a = {a}, b = {b}, x = {x}");
                // The reflected form I(1-x)(b, a) = 1 - Ix(a, b).
                Assert.AreEqual(1d - expected, Beta.Incomplete(b, a, 1d - x), 1E-13, $"a = {a}, b = {b}, x = {x}");
            }
        }

        /// <summary>
        /// One shape from 30 to 1E4 and the other from 0.01 to 1, with x = 1 - c/a for c from 0.01 to 30, where the
        /// asymptotic expansion is used: Ix(a, b) matches mpmath (the incomplete beta function at 50 digits) on these 92
        /// points to 2.6E-14 relative, and I(1-x)(b, a) = 1 - Ix(a, b) to 2.3E-16 absolute.
        /// </summary>
        [TestMethod]
        public void Test_Incomplete_LargeShapeGrid_MpmathOracle()
        {
            var expected = new double[]
            {
                0.039808736918165652, 0.018268051888993071, 0.0022149011592203284, 0.00011580620598255764,
                5.0597185913837583E-9, 0.33835640664330524, 0.17346008181410169, 0.024054726773150899,
                0.0013589913283754034, 6.5313235464814441E-8, 0.88799230693479569, 0.65579482777073417,
                0.15551940585192561, 0.012282448499852760, 8.9912739638672047E-7, 0.99004818330071176,
                0.90468628845713883, 0.36166151346161057, 0.042391158275216237, 5.2150950508465805E-6,
                0.039699246333342729, 0.018174963599108463, 0.0022159865924809225, 0.00012786355560848262,
                2.5262475386365702E-8, 1.0972085851343874E-19, 0.33767126972025158, 0.17275084569668611,
                0.024107288705198120, 0.0015037563797772467, 3.2696922058138912E-7, 1.5599469136124848E-18,
                0.88767385107756203, 0.65504351632442959, 0.15677586542444077, 0.013697876866621616,
                4.5452549045170283E-6, 3.2903699569020347E-17, 0.99004933869137191, 0.90479214711370892,
                0.36603234127322920, 0.047552507925405628, 2.6561398887587544E-5, 3.2344765096247375E-16,
                0.039657236583211623, 0.018139274120032151, 0.0022162154994438497, 0.00013246237642282473,
                4.0749421054770045E-8, 1.9887494991380517E-17, 0.33740796311685561, 0.17247846418621696,
                0.024125403119185585, 0.0015591367720306690, 5.2794840861529427E-7, 2.8306269488103688E-16,
                0.88755076825435153, 0.65475313089791687, 0.15724727426672383, 0.014245064052179493,
                7.3663977380150530E-6, 5.9950522734167458E-15, 0.99004978424639267, 0.90483289355855623,
                0.36769542477096373, 0.049563082823155427, 4.3171247410657862E-5, 5.9119978336199198E-14,
                0.039653042419865926, 0.018135711791061420, 0.0022162327670916882, 0.00013292047708366470,
                4.2606939230088363E-8, 3.0065365990931709E-17, 0.33738166205122028, 0.17245126242347034,
                0.024127150309053508, 0.0015646584217833958, 5.5206991156403150E-7, 4.2797357882056965E-16,
                0.88753845248248420, 0.65472407468418625, 0.15729401776335150, 0.014299796921443782,
                7.7058115692948758E-6, 9.0678429533711715E-15, 0.99004982879863090, 0.90483696561475935,
                0.36786104643297046, 0.049764664747471460, 4.5173345977048242E-5, 8.9450570868922940E-14,
            };
            int index = 0;
            foreach (double a in new[] { 30d, 100d, 1000d, 1E4 })
            {
                foreach (double b in new[] { 0.01, 0.1, 0.5, 1d })
                {
                    foreach (double c in new[] { 0.01, 0.1, 1d, 3d, 10d, 30d })
                    {
                        double x = 1d - c / a;
                        if (!(x > 0d)) continue;
                        double value = expected[index++];
                        Assert.AreEqual(value, Beta.Incomplete(a, b, x), 5E-14 * value, $"a = {a}, b = {b}, x = {x}");
                        Assert.AreEqual(1d - value, Beta.Incomplete(b, a, 1d - x), 1E-15, $"a = {a}, b = {b}, x = {x}");
                    }
                }
            }
            Assert.AreEqual(expected.Length, index);
        }

        /// <summary>
        /// Small tail targets keep their relative accuracy. Reference values are mpmath roots at 50 digits of
        /// ln Ix(a, b) = ln y (or of the upper tail for y &gt; 1/2). The points include (100, 0.5, 1E-12), (300, 0.5, 1E-20),
        /// (300, 2, 1E-20), (1000, 100, 1E-20) and (10000, 20, 1E-300). The largest error measured on these points is
        /// 2.2E-16 relative.
        /// </summary>
        [TestMethod]
        public void Test_IncompleteInverse_SmallTail_MpmathOracle()
        {
            var cases = new (double a, double b, double y, double expected)[]
            {
                (100d, 0.5d, 1E-12, 0.77503783769810146),
                (100d, 0.5d, 2E-10, 0.81641606901002351),
                (300d, 0.5d, 1E-20, 0.86468569608042303),
                (1000d, 1d, 1E-100, 0.79432823472428150),
                (1E4, 0.1d, 1E-300, 0.93401025014305818),
                (300d, 2d, 1E-20, 0.84676040157092620),
                (1000d, 100d, 1E-20, 0.80863075850883868),
                (1E4, 20d, 1E-300, 0.92522316039532233),
                (0.1d, 100d, 0.999999999999, 0.20257619700514166),
                (100d, 1d, 1E-6, 0.87096358995608064),
                (75d, 0.5d, 0.999999, 0.99999999999998949),
            };
            foreach (var (a, b, y, expected) in cases)
            {
                Assert.AreEqual(expected, Beta.IncompleteInverse(a, b, y), 1E-15 * expected, $"a = {a}, b = {b}, y = {y}");
            }
        }

        /// <summary>
        /// A root below e^-708 comes from the leading series term (y a B(a, b))^(1/a): it is subnormal or 0 where the root
        /// is below the smallest positive double, where the inverse threw ArithmeticException before. mpmath roots:
        /// (0.5, 100, 1E-300) is 7.9E-603 and (0.1, 0.5, 1E-50) is 3.5E-500, both below the double range.
        /// </summary>
        [TestMethod]
        public void Test_IncompleteInverse_UnderflowingRoot()
        {
            Assert.AreEqual(0d, Beta.IncompleteInverse(0.5d, 100d, 1E-300), 0d);
            Assert.AreEqual(0d, Beta.IncompleteInverse(0.1d, 0.5d, 1E-50), 0d);
            Assert.AreEqual(0d, Beta.IncompleteInverse(0.1d, 0.5d, 1E-300), 0d);
            Assert.AreEqual(0d, Beta.IncompleteInverse(0.5d, 5d, 1E-300), 0d);
            // Subnormal roots: mpmath 4.4444444444444446E-311 and 1.3321008721895818E-311.
            Assert.AreEqual(4.4444444444444446E-311, Beta.IncompleteInverse(0.5d, 2d, 1E-155), 1E-12 * 4.4444444444444446E-311);
            Assert.AreEqual(1.3321008721895818E-311, Beta.IncompleteInverse(0.1d, 5d, 1E-31), 1E-12 * 1.3321008721895818E-311);
        }

        /// <summary>
        /// A small shape below 0.1 with the other shape between 1 and 15 and x close to 1, where the continued fractions
        /// converge slowly near x = a / (a + b): the large shape is raised to 15 or more by the upward recurrence before the
        /// asymptotic expansion. mpmath at 60 digits; the points include (14.999, 1E-10, 0.999999) and (10, 1E-10,
        /// 1 - 1E-12). The largest error measured on these points is 1.8E-15 relative.
        /// </summary>
        [TestMethod]
        public void Test_Incomplete_SmallShapeBelowAsymptoticShape_MpmathOracle()
        {
            var cases = new (double a, double b, double x, double expected)[]
            {
                (14.999, 1E-10, 0.999999, 1.0564031165429819E-9),
                (14.999, 1E-05, 0.999999, 0.00010563481075749502),
                (10d, 1E-05, 0.999999, 0.00010985955502365953),
                (5d, 1E-10, 0.999999, 1.1732181217788162E-9),
                (1.5, 1E-10, 0.999999, 1.3201805410376596E-9),
                (10d, 1E-10, 0.999999999999, 2.4802074953253948E-9),
                (3d, 1E-10, 0.999999999999, 2.6131043203816287E-9),
            };
            foreach (var (a, b, x, expected) in cases)
            {
                Assert.AreEqual(expected, Beta.Incomplete(a, b, x), 1E-14 * expected, $"a = {a}, b = {b}, x = {x}");
                Assert.AreEqual(1d - expected, Beta.Incomplete(b, a, 1d - x), 1E-15, $"a = {a}, b = {b}, x = {x}");
            }
        }

        /// <summary>
        /// A large shape above 2^996, where the double-double residual of ν ln x needs scaling, gives finite values:
        /// Ix(1E305, 1/2) underflows to +0 for x ≤ 1 - 1E-9, and Ix(1/2, 1E305) at x = 1E-305 is erf(1) = 0.84270079294971487
        /// (the limit P(1/2, b x)).
        /// </summary>
        [TestMethod]
        public void Test_Incomplete_HugeShape_Finite()
        {
            foreach (double x in new[] { 0.7d, 0.8d, 0.9d, 1d - 1E-9 })
            {
                double value = Beta.Incomplete(1E305, 0.5d, x);
                Assert.AreEqual(0d, value, 0d, $"x = {x}");
                Assert.AreEqual(0L, System.BitConverter.DoubleToInt64Bits(value), $"x = {x} returns +0");
                Assert.AreEqual(0d, Beta.Incomplete(1.7E308, 0.5d, x), 0d, $"x = {x}");
            }
            Assert.AreEqual(0.84270079294971487, Beta.Incomplete(0.5d, 1E305, 1E-305), 1E-15);
        }
    }
}
