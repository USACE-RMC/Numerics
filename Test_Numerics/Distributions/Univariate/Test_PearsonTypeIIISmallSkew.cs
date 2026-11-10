using System;
using Microsoft.VisualStudio.TestTools.UnitTesting;
using Numerics;
using Numerics.Distributions;

namespace Distributions.Univariate
{
    /// <summary>
    /// Testing the Pearson Type III conditional moments and the small-skew density and tails.
    /// </summary>
    [TestClass]
    public class Test_PearsonTypeIIISmallSkew
    {
        /// <summary>
        /// Standardized conditional moments of PIII(0, 1, γ): γ, lower limit, upper limit, E[Z], E[Z²], E[Z³], E[Z⁴].
        /// </summary>
        /// <remarks>
        /// mpmath 1.4.1 at 60 digits: mpmath.quad(z^k f(z), [a, b]) / mpmath.quad(f(z), [a, b]) with the standardized density
        /// log f(z) = (α-1) log(α + 2z/γ) - (α + 2z/γ) - loggamma(α) + log(2/γ), α = 4/γ², for γ &gt; 0, f(z; γ) = f(-z; -γ)
        /// for γ &lt; 0, and the standard Normal density at γ = 0. Infinite limits are truncated at |z| = 60 and at the support.
        /// </remarks>
        private static readonly double[,] QuadratureMoments =
        {
            { 0.0d, -1.0d, 1.0d, 0.0d, 0.2911250947727932d, 0.0d, 0.16450037909117285d },
            { 0.0d, -1.0d, double.PositiveInfinity, 0.2875999709391784d, 0.7124000290608217d, 0.8627999128175351d, 1.8496001162432865d },
            { 0.0d, double.NegativeInfinity, 0.5d, -0.5091604338370335d, 0.7454197830814833d, -1.1456109761333253d, 2.1726142950148204d },
            { 0.0d, 2.0d, 5.0d, 2.3731800849531375d, 5.746164117648071d, 14.237708143911394d, 36.216286994499335d },
            { 1e-06d, -1.0d, 1.0d, -1.1814581753787273e-07d, 0.2911250947728043d, -6.331235784081389e-08d, 0.16450037909118131d },
            { 1e-06d, -1.0d, double.PositiveInfinity, 0.28759992300582643d, 0.712400220794135d, 0.8628004814177002d, 1.849602033577301d },
            { 1e-06d, double.NegativeInfinity, 0.5d, -0.5091604120389883d, 0.7454195394002998d, -1.1456101816681843d, 2.1726118482807535d },
            { 1e-06d, 2.0d, 5.0d, 2.3731804324698373d, 5.746165995419823d, 14.237715948214545d, 36.21631661429944d },
            { -1e-06d, -1.0d, 1.0d, 1.1814581753787273e-07d, 0.2911250947728043d, 6.331235784081389e-08d, 0.16450037909118131d },
            { -1e-06d, -1.0d, double.PositiveInfinity, 0.28760001887248343d, 0.7123998373275071d, 0.862799344217613d, 1.8495981989110217d },
            { -1e-06d, double.NegativeInfinity, 0.5d, -0.5091604556350571d, 0.7454200267626993d, -1.1456117705989053d, 2.1726167417513715d },
            { -1e-06d, 2.0d, 5.0d, 2.373179737436382d, 5.746162239876494d, 14.237700339611242d, 36.216257374720364d },
            { 0.0001d, -1.0d, 1.0d, -1.181458175872053e-05d, 0.29112509488371485d, -6.331235787790104e-06d, 0.16450037917590324d },
            { 0.0001d, -1.0d, double.PositiveInfinity, 0.287595177372139d, 0.7124192023867295d, 0.8628567740366558d, 1.849791858304155d },
            { 0.0001d, double.NegativeInfinity, 0.5d, -0.5091582539259856d, 0.745395415124311d, -1.1455315317919552d, 2.1723696339024157d },
            { 0.0001d, 2.0d, 5.0d, 2.3732148363482755d, 5.746351895692999d, 14.238488589061632d, 36.21924907908413d },
            { -0.0001d, -1.0d, 1.0d, 1.181458175872053e-05d, 0.29112509488371485d, 6.331235787790104e-06d, 0.16450037917590324d },
            { -0.0001d, -1.0d, double.PositiveInfinity, 0.28760476403784513d, 0.7123808557239529d, 0.862743054027963d, 1.8494083916759094d },
            { -0.0001d, double.NegativeInfinity, 0.5d, -0.5091626135328609d, 0.7454441513642462d, -1.1456904248640734d, 2.1728589809648606d },
            { -0.0001d, 2.0d, 5.0d, 2.373145333002781d, 5.745976341360421d, 14.23692772873311d, 36.21332512118304d },
            { 0.001d, -1.0d, 1.0d, -0.00011814582247161893d, 0.2911251058649607d, -6.331236154989607e-05d, 0.16450038756421176d },
            { 0.001d, -1.0d, double.PositiveInfinity, 0.28755201418799275d, 0.7125917618191012d, 0.8633686343257975d, 1.8515183242207995d },
            { 0.001d, double.NegativeInfinity, 0.5d, -0.5091386250424771d, 0.7451761181662402d, -1.1448167302274073d, 2.17016880127307d },
            { 0.001d, 2.0d, 5.0d, 2.3735275739020296d, 5.748041977010489d, 14.245513943325541d, 36.24591734452237d },
            { -0.001d, -1.0d, 1.0d, 0.00011814582247161893d, 0.2911251058649607d, 6.331236154989607e-05d, 0.16450038756421176d },
            { -0.001d, -1.0d, double.PositiveInfinity, 0.28764788085310944d, 0.712208295206464d, 0.8622314342641217d, 1.8476836576148865d },
            { -0.001d, double.NegativeInfinity, 0.5d, -0.5091822211095353d, 0.7456634805557871d, -1.1464056609770101d, 2.175062272520135d },
            { -0.001d, 2.0d, 5.0d, 2.3728325404823205d, 5.74428643401291d, 14.229905341689483d, 36.18667777130666d },
            { 0.01d, -1.0d, 1.0d, -0.001181463109175591d, 0.2911262039993448d, -0.0006331272875333255d, 0.16450122640361198d },
            { 0.01d, -1.0d, double.PositiveInfinity, 0.2871182917122436d, 0.7143172997463176d, 0.868498048134194d, 1.8688610782487223d },
            { 0.01d, double.NegativeInfinity, 0.5d, -0.5089413782457515d, 0.7429846039858955d, -1.1376882550130818d, 2.148270815851771d },
            { 0.01d, 2.0d, 5.0d, 2.3766524580173765d, 5.764950452350169d, 14.315900170805065d, 36.51353821073943d },
            { -0.01d, -1.0d, 1.0d, 0.001181463109175591d, 0.2911262039993448d, 0.0006331272875333255d, 0.16450122640361198d },
            { -0.01d, -1.0d, double.PositiveInfinity, 0.2880769664190304d, 0.7104826487488745d, 0.8571260727696025d, 1.830514088736049d },
            { -0.01d, double.NegativeInfinity, 0.5d, -0.509377337221865d, 0.7478582180751768d, -1.153577590929948d, 2.1972061509367466d },
            { -0.01d, 2.0d, 5.0d, 2.3697021590535083d, 5.727395350557649d, 14.159815803536778d, 35.921148273638096d },
            { 0.02d, -1.0d, 1.0d, -0.002362955822242721d, 0.2911295317978147d, -0.0012662768307718677d, 0.16450376844455813d },
            { 0.02d, -1.0d, double.PositiveInfinity, 0.2866319041952485d, 0.716234414846704d, 0.8742204008826795d, 1.888297952371344d },
            { 0.02d, double.NegativeInfinity, 0.5d, -0.5087201755769495d, 0.7405527104557558d, -1.1298093408390213d, 2.124173829194978d },
            { 0.02d, 2.0d, 5.0d, 2.3801191676185747d, 5.783753329389021d, 14.39438668962903d, 36.812883226356284d },
            { -0.02d, -1.0d, 1.0d, 0.002362955822242721d, 0.2911295317978147d, 0.0012662768307718677d, 0.16450376844455813d },
            { -0.02d, -1.0d, double.PositiveInfinity, 0.2885493024325924d, 0.7085652045430817d, 0.8514766032069156d, 1.811602013100445d },
            { -0.02d, double.NegativeInfinity, 0.5d, -0.509592083259364d, 0.7502998792029116d, -1.1615881849176273d, 2.222048272748843d },
            { -0.02d, 2.0d, 5.0d, 2.366218783160299d, 5.708645114415118d, 14.082227948156353d, 35.62813847120185d },
            { 0.59d, -1.0d, 1.0d, -0.07075388941215192d, 0.29511153612034924d, -0.03814586192544972d, 0.16755945405395878d },
            { 0.59d, -1.0d, double.PositiveInfinity, 0.2501660121202392d, 0.8236329614552313d, 1.236441483619304d, 3.314983585248539d },
            { 0.59d, double.NegativeInfinity, 0.5d, -0.4927283233662496d, 0.6082809829238316d, -0.749752947649001d, 1.0997205496813476d },
            { 0.59d, 2.0d, 5.0d, 2.563053993380622d, 6.834309428210231d, 19.075270116739695d, 56.02089430765677d },
            { -0.59d, -1.0d, 1.0d, 0.07075388941215192d, 0.29511153612034924d, 0.03814586192544972d, 0.16755945405395878d },
            { -0.59d, -1.0d, double.PositiveInfinity, 0.30845536232992904d, 0.6005503057827418d, 0.5710414065779694d, 0.9878239101967936d },
            { -0.59d, double.NegativeInfinity, 0.5d, -0.5180891510618184d, 0.8937917240323272d, -1.6930377070681644d, 4.11495239896958d },
            { -0.59d, 2.0d, 5.0d, 2.1594038754992373d, 4.681783607726199d, 10.194170924436966d, 22.298740559045783d },
            { 0.6d, -1.0d, 1.0d, -0.07199083845507898d, 0.29525244713684695d, -0.03882104708312876d, 0.16766809770909552d },
            { 0.6d, -1.0d, double.PositiveInfinity, 0.24935251526467403d, 0.8254532393147281d, 1.243329489382859d, 3.3460037431240837d },
            { 0.6d, double.NegativeInfinity, 0.5d, -0.49239139362253115d, 0.6060868851019751d, -0.7442285045895101d, 1.0869060769725498d },
            { 0.6d, 2.0d, 5.0d, 2.5659179941564068d, 6.8515166429386465d, 19.155754745354006d, 56.36837815166333d },
            { -0.6d, -1.0d, 1.0d, 0.07199083845507898d, 0.29525244713684695d, 0.03882104708312876d, 0.16766809770909552d },
            { -0.6d, -1.0d, double.PositiveInfinity, 0.308692028571704d, 0.5987003628567847d, 0.5668558680010413d, 0.9772387787977131d },
            { -0.6d, double.NegativeInfinity, 0.5d, -0.5181691814725278d, 0.8963661637054945d, -1.703700356536484d, 4.157657664315253d },
            { -0.6d, 2.0d, 5.0d, 2.155625793746465d, 4.6645638493689905d, 10.135016452857396d, 22.117183090507037d },
        };

        /// <summary>
        /// Small-skew log density and log tails of PIII(0, 1, γ): γ, x, log f(x), log F(x), log(1 - F(x)).
        /// </summary>
        /// <remarks>
        /// R 4.4.3 with exactly representable gamma coordinates a = 4/γ² and y = a + 2x/γ: dgamma(y, a, log = TRUE) + log(2/|γ|),
        /// pgamma(y, a, log.p = TRUE) and pgamma(y, a, lower.tail = FALSE, log.p = TRUE), with the two tails exchanged for γ &lt; 0.
        /// </remarks>
        private static readonly double[,] SmallSkewTails =
        {
            { 1e-06d, -30.0d, -450.9234235837178d, -454.3257439903756d, -4.884683154279846e-198d },
            { 1e-06d, -20.0d, -200.92026187648813d, -203.91848868967645d, -2.749955106933085e-89d },
            { 1e-06d, -5.0d, -13.418956866573966d, -15.065019140043667d, -2.8664566613448664e-07d },
            { 1e-06d, 0.5d, -1.043938762371333d, -0.3689463516436048d, -1.1759119042283486d },
            { 1e-06d, 3.0d, -5.41893553320863d, -0.001350815881873037d, -6.607721844050543d },
            { 1e-06d, 5.0d, -13.4189201999073d, -2.866575598940146e-07d, -15.06497764801193d },
            { 1e-06d, 20.0d, -200.9176152098213d, -2.757297970357437e-89d, -203.91582207251813d },
            { 1e-06d, 30.0d, -450.9144535837166d, -4.928843563598553e-198d, -454.31674402356083d },
            { -1e-06d, -30.0d, -450.9144535837166d, -454.31674402356083d, -4.928843563598553e-198d },
            { -1e-06d, -20.0d, -200.9176152098213d, -203.91582207251813d, -2.757297970357437e-89d },
            { -1e-06d, -5.0d, -13.4189201999073d, -15.06497764801193d, -2.866575598940146e-07d },
            { -1e-06d, 0.5d, -1.043938304037999d, -0.36894647893371324d, -1.175911618958906d },
            { -1e-06d, 3.0d, -5.418941533208633d, -0.0013508040476355896d, -6.60773059898029d },
            { -1e-06d, 5.0d, -13.418956866573966d, -2.8664566613448664e-07d, -15.065019140043667d },
            { -1e-06d, 20.0d, -200.92026187648813d, -2.749955106933085e-89d, -203.91848868967645d },
            { -1e-06d, 30.0d, -450.9234235837178d, -4.884683154279846e-198d, -454.3257439903756d },
            { 0.0001d, -30.0d, -451.36794426554707d, -454.7717491554914d, -3.1270787352673867e-198d },
            { 0.0001d, -20.0d, -201.05137144647946d, -204.05058630926894d, -2.409663153016309e-89d },
            { 0.0001d, -5.0d, -13.420772226194266d, -15.067073386389279d, -2.860574296363213e-07d },
            { 0.0001d, 0.5d, -1.0439614498062397d, -0.36894005081011655d, -1.1759260251526142d },
            { 0.0001d, 3.0d, -5.418638572783056d, -0.0013514017378129864d, -6.607288525685341d },
            { 0.0001d, 5.0d, -13.417105559381774d, -2.87246806022273e-07d, -15.062924183059343d },
            { 0.0001d, 20.0d, -200.7867046204794d, -3.1460536031891e-89d, -203.78392443345749d },
            { 0.0001d, 30.0d, -450.470943052795d, -7.69135640630715e-198d, -453.871751259122d },
            { -0.0001d, -30.0d, -450.470943052795d, -453.871751259122d, -7.69135640630715e-198d },
            { -0.0001d, -20.0d, -200.7867046204794d, -203.78392443345749d, -3.1460536031891e-89d },
            { -0.0001d, -5.0d, -13.417105559381774d, -15.062924183059343d, -2.87246806022273e-07d },
            { -0.0001d, 0.5d, -1.0439156164728978d, -0.3689527798209626d, -1.1758974982083539d },
            { -0.0001d, 3.0d, -5.419238572792956d, -0.0013502183140740999d, -6.608164018672132d },
            { -0.0001d, 5.0d, -13.420772226194266d, -2.860574296363213e-07d, -15.067073386389279d },
            { -0.0001d, 20.0d, -201.05137144647946d, -2.409663153016309e-89d, -204.05058630926894d },
            { -0.0001d, 30.0d, -451.36794426554707d, -3.1270787352673867e-198d, -454.7717491554914d },
            { 0.001d, -30.0d, -455.4550651089107d, -458.87246757628577d, -5.178684272311935e-200d },
            { 0.001d, -20.0d, -202.25230222394924d, -205.26054464251177d, -7.185849370522326e-90d },
            { 0.001d, -5.0d, -13.437307897941333d, -15.08578356170015d, -2.807550035725487e-07d },
            { 0.001d, 0.5d, -1.044167693365349d, -0.3688827729227111d, -1.176054405002007d },
            { 0.001d, 3.0d, -5.415942486594325d, -0.0013567332001637369d, -6.603353817411734d },
            { 0.001d, 5.0d, -13.40064108544067d, -2.926491963308346e-07d, -15.04429137371614d },
            { 0.001d, 20.0d, -199.60547621255978d, -1.0342956984650212e-88d, -202.59376747296827d },
            { 0.001d, 30.0d, -446.48385216391245d, -4.20129436319156e-196d, -449.8712855673096d },
            { -0.001d, -30.0d, -446.48385216391245d, -449.8712855673096d, -4.20129436319156e-196d },
            { -0.001d, -20.0d, -199.60547621255978d, -202.59376747296827d, -1.0342956984650212e-88d },
            { -0.001d, -5.0d, -13.40064108544067d, -15.04429137371614d, -2.926491963308346e-07d },
            { -0.001d, 0.5d, -1.0437093600231613d, -0.3690100630312332d, -1.1757691355583377d },
            { -0.001d, 3.0d, -5.421942496494342d, -0.0013448989683984647d, -6.612108759301565d },
            { -0.001d, 5.0d, -13.437307897941333d, -2.807550035725487e-07d, -15.08578356170015d },
            { -0.001d, 20.0d, -202.25230222394924d, -7.185849370522326e-90d, -205.26054464251177d },
            { -0.001d, 30.0d, -455.4550651089107d, -5.178684272311935e-200d, -458.87246757628577d },
        };

        /// <summary>
        /// Conditional moments of PIII(100, 15, γ) for |γ| at or above 0.6: γ, a, b, E[X], E[(X-µ)²], E[(X-µ)³], E[(X-µ)⁴].
        /// </summary>
        /// <remarks>Values of the truncated-gamma moment path, pinned bit for bit.</remarks>
        private static readonly double[,] TruncatedGammaMoments =
        {
            { 0.6d, 85d, 115d, 98.9201374233989d, 66.4318006008791d, -131.02103386988077d, 8488.197445749547d },
            { 0.6d, 85d, double.PositiveInfinity, 103.74028772855824d, 185.72697883918514d, 4196.23702726927d, 169391.43956188593d },
            { 0.6d, double.NegativeInfinity, 107.5d, 92.61412909632557d, 136.36954914505725d, -2511.771203222097d, 55024.62016176899d },
            { 0.6d, 130d, 175d, 138.4887699149619d, 1541.5912449068942d, 64650.67228352937d, 2853649.1451394456d },
            { 0.6d, 104.5d, 110.5d, 107.36184266878762d, 57.17429524639013d, 465.2408406449171d, 3936.1217247059026d },
            { 0.6d, double.NegativeInfinity, 82d, 77.24515605373259d, 532.1903933505153d, -12822.535677993928d, 318830.5338957446d },
            { 0.6d, 122.5d, double.PositiveInfinity, 131.76479679315042d, 1082.6495135490743d, 40118.93256592272d, 1634214.900439044d },
            { 0.6d, 62.5d, 92.5d, 84.58526134912299d, 270.7918840249892d, -5346.229549224721d, 116324.76187877089d },
            { -0.6d, 85d, 115d, 101.07986257660109d, 66.4318006008791d, 131.02103386988077d, 8488.197445746562d },
            { -0.6d, 85d, double.PositiveInfinity, 104.63038042802592d, 134.70758164004192d, 1913.138554706451d, 49472.71318936557d },
            { -0.6d, double.NegativeInfinity, 107.5d, 92.22746227879082d, 201.6823868249227d, -5749.988703939837d, 210481.4193371861d },
            { -0.6d, 130d, 175d, 132.33438691017386d, 1049.526866380939d, 34205.68054251386d, 1119682.3946097787d },
            { -0.6d, 104.5d, 110.5d, 107.4525582918026d, 58.51914759537185d, 480.6804223182025d, 4098.381785721182d },
            { -0.6d, double.NegativeInfinity, 82d, 72.2545964460988d, 849.2715800565394d, -29118.386578372367d, 1128167.7293268216d },
            { -0.6d, 122.5d, double.PositiveInfinity, 126.22242374005467d, 697.0036274291364d, 18802.160060443955d, 515338.0833563858d },
            { -0.6d, 62.5d, 92.5d, 82.80392763621003d, 350.7596930887321d, -8235.6559851539d, 214664.3923144173d },
            { 0.65d, 85d, 115d, 98.8269044336371d, 66.59951182387258d, -142.4942670305465d, 8517.32118640731d },
            { 0.65d, 85d, double.PositiveInfinity, 103.6777184121384d, 187.76310107471204d, 4313.1501638821d, 177407.6147456083d },
            { 0.65d, double.NegativeInfinity, 107.5d, 92.63980985748766d, 133.91764698401704d, -2420.3992016691755d, 51890.993177635704d },
            { 0.65d, 130d, 175d, 138.70040985742796d, 1560.7263013546499d, 65997.83304042669d, 2941204.056634417d },
            { 0.65d, 104.5d, 110.5d, 107.3590479883259d, 57.133104903359005d, 464.7705395115766d, 3931.204033350598d },
            { 0.65d, double.NegativeInfinity, 82d, 77.46390329676285d, 520.7862692465231d, -12365.272723952547d, 302119.13415680296d },
            { 0.65d, 122.5d, double.PositiveInfinity, 131.98149966868223d, 1100.4935534599792d, 41312.12120678781d, 1711312.1910558415d },
            { 0.65d, 62.5d, 92.5d, 84.72107209415842d, 264.95574294143717d, -5144.82383217514d, 109783.79724890941d },
            { -0.75d, 85d, 115d, 101.36213319399162d, 66.98267072006166d, 165.88486034472083d, 8584.0464534332d },
            { -0.75d, 85d, double.PositiveInfinity, 104.67734936901046d, 128.52966926401982d, 1711.2520448671974d, 42094.094375809305d },
            { -0.75d, double.NegativeInfinity, 107.5d, 92.21403207118773d, 210.40131013335503d, -6308.6610029744115d, 245194.8335462139d },
            { -0.75d, 130d, 175d, 131.4744775092958d, 992.1903892906539d, 31328.40275812848d, 990872.6089837217d },
            { -0.75d, 104.5d, 110.5d, 107.4684914144776d, 58.75692995767467d, 483.42736761868275d, 4127.4163655071925d },
            { -0.75d, double.NegativeInfinity, 82d, 71.66518747109633d, 894.4099459963538d, -31993.256789909374d, 1308861.5485514211d },
            { -0.75d, 122.5d, double.PositiveInfinity, 125.47943656073119d, 654.9654919633996d, 16996.34942666603d, 445515.0175770388d },
            { -0.75d, 62.5d, 92.5d, 82.70229872728167d, 355.45465579681684d, -8409.830843692082d, 220737.93787053722d },
            { 0.8d, 85d, 115d, 98.54188676888656d, 67.19939766413727d, -177.8336590317117d, 8621.904897138738d },
            { 0.8d, 85d, double.PositiveInfinity, 103.47351979086055d, 193.73832188222528d, 4669.485721420423d, 203100.98096204756d },
            { 0.8d, double.NegativeInfinity, 107.5d, 92.72086338603937d, 126.73165571150282d, -2164.283042280259d, 43515.88708522195d },
            { 0.8d, 130d, 175d, 139.30338839418303d, 1615.774159003248d, 69913.34280551241d, 3198382.065999341d },
            { 0.8d, 104.5d, 110.5d, 107.35099889962883d, 57.014550536150594d, 463.4177703352355d, 3917.067039183359d },
            { 0.8d, double.NegativeInfinity, 82d, 78.12463107824175d, 487.50442706125705d, -11081.48242071302d, 257175.956245459d },
            { 0.8d, 122.5d, double.PositiveInfinity, 132.6256496599958d, 1154.831015310167d, 45056.249671119636d, 1962149.9700746571d },
            { 0.8d, 62.5d, 92.5d, 85.14367459999217d, 247.28448810005614d, -4553.600876555322d, 91219.7259677348d },
            { 1d, 85d, 115d, 98.14586122445964d, 68.25295286175995d, -227.74938056325945d, 8807.155696628544d },
            { 1d, 85d, double.PositiveInfinity, 103.15789473684208d, 201.3157894736833d, 5151.315789473807d, 241134.86842104053d },
            { 1d, double.NegativeInfinity, 107.5d, 92.83781400923098d, 117.56721013846476d, -1862.348505749821d, 34433.47824923341d },
            { 1d, 130d, 175d, 140.03350624069876d, 1683.3782263941787d, 74793.68976949403d, 3523786.280407466d },
            { 1d, 104.5d, 110.5d, 107.34098657161981d, 56.867244998740006d, 461.73870237525404d, 3899.537020177779d },
            { 1d, double.NegativeInfinity, 82d, 79.01455279359149d, 445.3471956672898d, -9562.528202750855d, 207839.60062130095d },
            { 1d, 122.5d, double.PositiveInfinity, 133.471189591078d, 1229.1356877323392d, 50443.860362453896d, 2345911.216310343d },
            { 1d, 62.5d, 92.5d, 85.71032779906123d, 225.00000000000006d, -3859.1465517252354d, 71072.65804595161d },
            { -1d, 85d, 115d, 101.85413877554036d, 68.25295286175995d, 227.74938056325945d, 8807.155696628184d },
            { -1d, 85d, double.PositiveInfinity, 104.73091057406762d, 118.55451208347795d, 1415.0469562434253d, 32218.915953387444d },
            { -1d, double.NegativeInfinity, 107.5d, 92.21153846153848d, 224.99999999999903d, -7317.908653846287d, 313242.1874999862d },
            { -1d, 130d, 175d, double.NaN, double.NaN, double.NaN, double.NaN },
            { -1d, 104.5d, 110.5d, 107.49858148049859d, 59.20728188097037d, 488.6441078269058d, 4182.694009778948d },
            { -1d, double.NegativeInfinity, 82d, 70.69816312083167d, 972.1968404187896d, -37262.574350758696d, 1665529.1028540155d },
            { -1d, 122.5d, double.PositiveInfinity, 124.2176902404785d, 588.2653536071808d, 14334.185988349926d, 350414.5568424267d },
            { -1d, 62.5d, 92.5d, 82.55981790110599d, 362.06002931954305d, -8655.571029984514d, 229327.1078900525d },
            { 1.5d, 85d, 115d, 97.01703695123985d, 72.70387113868864d, -377.66295729263203d, 9612.984433368047d },
            { 1.5d, 85d, double.PositiveInfinity, 102.06825948363779d, 217.24402693635827d, 6284.065757523561d, 351746.56174118526d },
            { 1.5d, double.NegativeInfinity, 107.5d, 93.16836172557966d, 96.9067823546187d, -1278.1142734463806d, 19393.623963531234d },
            { 1.5d, 130d, 175d, 141.52599728991981d, 1824.3006290447531d, 85173.95078263244d, 4229975.378678904d },
            { 1.5d, 104.5d, 110.5d, 107.31900697875032d, 56.54451292797512d, 458.0668678268249d, 3861.2672244733144d },
            { 1.5d, double.NegativeInfinity, 82d, 81.258094823834d, 351.50785993912046d, -6597.307757722269d, 123911.45962317978d },
            { 1.5d, 122.5d, double.PositiveInfinity, 135.5240925289588d, 1423.9381228523653d, 66008.52124499496d, 3593587.44140658d },
            { 1.5d, 62.5d, 92.5d, 87.10735172368213d, 176.65256896380782d, -2552.2203882002395d, 38542.13194038366d },
            { -2d, 85d, 115d, 104.69552928248997d, 84.13412152530097d, 645.4586199217043d, 11897.482804697975d },
            { -2d, 85d, double.PositiveInfinity, 104.69552928248997d, 84.13412152530097d, 645.4586199217043d, 11897.482804697975d },
            { -2d, double.NegativeInfinity, 107.5d, 92.5d, 281.25d, -12234.375000000005d, 737226.5624999987d },
            { -2d, 130d, 175d, double.NaN, double.NaN, double.NaN, double.NaN },
            { -2d, 104.5d, 110.5d, 107.69946869031841d, 62.257969645223596d, 524.4719256690257d, 4567.206921156069d },
            { -2d, double.NegativeInfinity, 82d, 67.0d, 1313.9999999999989d, -64962.00000000001d, 4002695.9999999893d },
            { -2d, 122.5d, double.PositiveInfinity, double.NaN, double.NaN, double.NaN, double.NaN },
            { -2d, 62.5d, 92.5d, 82.19552928248996d, 379.0853038132522d, -9292.88448525446d, 251712.61071139717d },
            { 3d, 85d, 115d, 95.58846771364604d, 59.56753926172646d, -297.2450265112364d, 5255.128245712297d },
            { 3d, 85d, double.PositiveInfinity, 100.0d, 225.00000000000006d, 10125.000000000007d, 835312.5000000005d },
            { 3d, double.NegativeInfinity, 107.5d, 94.30635748148035d, 54.190724444410336d, -443.82392500211415d, 4218.618624833834d },
            { 3d, 130d, 175d, 144.1545550427853d, 2079.481331712559d, 104503.86439128962d, 5581455.84097772d },
            { 3d, 104.5d, 110.5d, 107.27115688990412d, 55.84494774305825d, 450.1398516497143d, 3778.953641968424d },
            { 3d, double.NegativeInfinity, 82d, double.NaN, double.NaN, double.NaN, double.NaN },
            { 3d, 122.5d, double.PositiveInfinity, 141.2609968699066d, 2081.744859145807d, 133134.34691840954d, 10861734.73938734d },
            { 3d, 62.5d, 92.5d, 90.74532736651656d, 86.17991049774847d, -807.0820483023068d, 7597.71634282532d },
        };

        /// <summary>
        /// Verifies the standardized conditional moments against high-precision quadrature for small and moderate skew,
        /// including zero skew, where the shape is infinite.
        /// </summary>
        [TestMethod]
        public void Test_ConditionalMoments_SmallSkew_MatchesQuadrature()
        {
            for (int i = 0; i < QuadratureMoments.GetLength(0); i++)
            {
                double g = QuadratureMoments[i, 0];
                var actual = new PearsonTypeIII(0d, 1d, g).ConditionalMoments(QuadratureMoments[i, 1], QuadratureMoments[i, 2]);
                // At |γ| = 0.6 the truncated-gamma path applies; it agrees with quadrature to about 1E-9.
                double tolerance = Math.Abs(g) < 0.6d ? 2E-12 : 5E-9;
                for (int k = 0; k < 4; k++)
                {
                    double expected = QuadratureMoments[i, 3 + k];
                    Assert.AreEqual(expected, actual[k], tolerance * Math.Abs(expected) + 1E-15, $"g = {g}, row {i}, moment {k + 1}");
                }
            }
        }

        /// <summary>
        /// Verifies location and scale handling: the physical moments are µ + σE[Z] and σᵏE[Zᵏ].
        /// </summary>
        [TestMethod]
        public void Test_ConditionalMoments_LocationScale()
        {
            for (int i = 0; i < QuadratureMoments.GetLength(0); i++)
            {
                double g = QuadratureMoments[i, 0];
                if (Math.Abs(g) >= 0.6d) continue;
                double lower = QuadratureMoments[i, 1], upper = QuadratureMoments[i, 2];
                double a = double.IsInfinity(lower) ? lower : 100d + 15d * lower;
                double b = double.IsInfinity(upper) ? upper : 100d + 15d * upper;
                var actual = new PearsonTypeIII(100d, 15d, g).ConditionalMoments(a, b);
                Assert.AreEqual(100d + 15d * QuadratureMoments[i, 3], actual[0], 1E-12 * 100d);
                double scale = 15d;
                for (int k = 1; k < 4; k++)
                {
                    scale *= 15d;
                    double expected = scale * QuadratureMoments[i, 3 + k];
                    Assert.AreEqual(expected, actual[k], 2E-12 * Math.Abs(expected) + 1E-15 * scale);
                }
            }
        }

        /// <summary>
        /// Verifies the whole-line limit, which reproduces the unconditional moments µ, σ², γσ³ and (3 + 3γ²/2)σ⁴.
        /// </summary>
        [TestMethod]
        public void Test_ConditionalMoments_WholeLine()
        {
            foreach (double g in new[] { -0.5d, -0.01d, -1E-7d, 0d, 1E-7d, 0.01d, 0.5d })
            {
                var moments = new PearsonTypeIII(5d, 2d, g).ConditionalMoments(double.NegativeInfinity, double.PositiveInfinity);
                Assert.AreEqual(5d, moments[0], 1E-14);
                Assert.AreEqual(4d, moments[1], 1E-14);
                Assert.AreEqual(8d * g, moments[2], 1E-13);
                Assert.AreEqual(16d * (3d + 1.5d * g * g), moments[3], 1E-13);
            }
        }

        /// <summary>
        /// Verifies that zero skew reproduces the truncated Normal moments and that small skews approach them continuously.
        /// </summary>
        [TestMethod]
        public void Test_ConditionalMoments_ZeroSkewMatchesNormal()
        {
            var normal = new Normal(3d, 0.5d);
            foreach (var interval in new[] { new[] { 2d, 3.2d }, new[] { double.NegativeInfinity, 2.4d }, new[] { 3.1d, double.PositiveInfinity }, new[] { 3.5d, 4.5d } })
            {
                var expected = normal.ConditionalMoments(interval[0], interval[1]);
                var zero = new PearsonTypeIII(3d, 0.5d, 0d).ConditionalMoments(interval[0], interval[1]);
                var tiny = new PearsonTypeIII(3d, 0.5d, 1E-12).ConditionalMoments(interval[0], interval[1]);
                for (int k = 0; k < 4; k++)
                {
                    Assert.AreEqual(expected[k], zero[k], 1E-14 * Math.Max(1d, Math.Abs(expected[k])));
                    Assert.AreEqual(expected[k], tiny[k], 1E-11 * Math.Max(1d, Math.Abs(expected[k])));
                }
            }
        }

        /// <summary>
        /// Verifies the not-a-number contract for reversed, empty, out-of-support and negligible-mass intervals.
        /// </summary>
        [TestMethod]
        public void Test_ConditionalMoments_UndefinedIntervals()
        {
            var positive = new PearsonTypeIII(0d, 1d, 0.5d);
            var negative = new PearsonTypeIII(0d, 1d, -0.5d);
            var cases = new[]
            {
                positive.ConditionalMoments(1d, 1d),
                positive.ConditionalMoments(2d, -2d),
                // The support of PIII(0, 1, 0.5) is z > -4.
                positive.ConditionalMoments(double.NegativeInfinity, -5d),
                // The support of PIII(0, 1, -0.5) is z < 4.
                negative.ConditionalMoments(4.5d, 6d),
                // The interval mass is far below 1E-15.
                new PearsonTypeIII(0d, 1d, 0.01d).ConditionalMoments(40d, 41d),
                positive.ConditionalMoments(double.NaN, 1d),
            };
            foreach (var moments in cases)
                foreach (double m in moments) Assert.IsTrue(double.IsNaN(m));
        }

        /// <summary>
        /// Verifies that the truncated-gamma conditional moments for absolute skew at or above 0.6 keep their values bit for bit.
        /// </summary>
        [TestMethod]
        public void Test_ConditionalMoments_LargerSkew_SameBits()
        {
            for (int i = 0; i < TruncatedGammaMoments.GetLength(0); i++)
            {
                var actual = new PearsonTypeIII(100d, 15d, TruncatedGammaMoments[i, 0]).ConditionalMoments(TruncatedGammaMoments[i, 1], TruncatedGammaMoments[i, 2]);
                for (int k = 0; k < 4; k++)
                    Assert.AreEqual(TruncatedGammaMoments[i, 3 + k], actual[k], 0d, $"row {i}, moment {k + 1}");
            }
        }

        /// <summary>
        /// Verifies small-skew log densities and log tails against R, including the far tails where the rounded unit gamma
        /// coordinate alone would perturb the result by several parts per billion.
        /// </summary>
        [TestMethod]
        public void Test_SmallSkew_DensityAndTails_MatchR()
        {
            for (int i = 0; i < SmallSkewTails.GetLength(0); i++)
            {
                var d = new PearsonTypeIII(0d, 1d, SmallSkewTails[i, 0]);
                double x = SmallSkewTails[i, 1];
                Assert.AreEqual(SmallSkewTails[i, 2], d.LogPDF(x), 1E-12 * Math.Max(1d, Math.Abs(SmallSkewTails[i, 2])), $"row {i} LogPDF");
                Assert.AreEqual(SmallSkewTails[i, 3], d.LogCDF(x), 1E-12 * Math.Max(1d, Math.Abs(SmallSkewTails[i, 3])), $"row {i} LogCDF");
                Assert.AreEqual(SmallSkewTails[i, 4], d.LogCCDF(x), 1E-12 * Math.Max(1d, Math.Abs(SmallSkewTails[i, 4])), $"row {i} LogCCDF");
            }

            // Regression at γ = 1E-6, x = 30 on the probability scale.
            // R: dgamma(4e12 + 6e7, 4e12, log = TRUE) + log(2e6) = -450.91445358371658
            // R: pgamma(4e12 + 6e7, 4e12, lower.tail = FALSE, log.p = TRUE) = -454.31674402356083
            var p3 = new PearsonTypeIII(0d, 1d, 1E-6);
            double pdf = Math.Exp(-450.91445358371658d), ccdf = Math.Exp(-454.31674402356083d);
            Assert.AreEqual(pdf, p3.PDF(30d), 1E-12 * pdf);
            Assert.AreEqual(ccdf, p3.CCDF(30d), 1E-12 * ccdf);
            // Location and scale.
            var scaled = new PearsonTypeIII(1000d, 50d, 1E-6);
            Assert.AreEqual(-450.91445358371658d - Math.Log(50d), scaled.LogPDF(2500d), 1E-12 * 450d);
            Assert.AreEqual(-454.31674402356083d, scaled.LogCCDF(2500d), 1E-12 * 454d);
        }

        /// <summary>
        /// Verifies the small-skew tails at the support endpoint, outside the support, at infinity and at not-a-number.
        /// </summary>
        [TestMethod]
        public void Test_SmallSkew_SupportLimits()
        {
            var d = new PearsonTypeIII(0d, 1d, 0.005d);
            double xi = d.Xi;
            Assert.AreEqual(double.NegativeInfinity, d.LogPDF(xi - 1d));
            Assert.AreEqual(double.NegativeInfinity, d.LogPDF(xi));
            Assert.AreEqual(0d, d.CDF(xi));
            Assert.AreEqual(1d, d.CCDF(xi));
            Assert.AreEqual(0d, d.PDF(double.PositiveInfinity));
            Assert.AreEqual(1d, d.CDF(double.PositiveInfinity));
            Assert.AreEqual(0d, d.CCDF(double.PositiveInfinity));
            Assert.IsTrue(double.IsNaN(d.LogPDF(double.NaN)));
            Assert.IsTrue(double.IsNaN(d.LogCDF(double.NaN)));
            var m = new PearsonTypeIII(0d, 1d, -0.005d);
            Assert.AreEqual(1d, m.CDF(m.Xi));
            Assert.AreEqual(0d, m.CDF(double.NegativeInfinity));
            // The two tails remain complementary across the small-skew region.
            foreach (double g in new[] { -0.009d, -1E-5d, 2E-6d, 0.004d })
                foreach (double x in new[] { -12.5d, -3.3d, 0.1d, 2.7d, 9.25d })
                {
                    var p = new PearsonTypeIII(0d, 1d, g);
                    Assert.AreEqual(1d, p.CDF(x) + p.CCDF(x), 1E-15);
                }
        }

        /// <summary>
        /// Verifies that the density, distribution and quantile functions keep their bits outside the small-skew region.
        /// </summary>
        [TestMethod]
        public void Test_DensityAndTails_OutsideSmallSkewRegion_SameBits()
        {
            double mu = 10d, sigma = 2d;
            foreach (double g in new[] { 0d, 0.1d, -0.1d, 0.1000000001d, 0.2d, -0.3d, 0.6d, 0.99d, -0.99d, 1d, 1.5d, -3d })
            {
                var d = new PearsonTypeIII(mu, sigma, g);
                foreach (double z in new[] { -40d, -8d, -3d, -1d, -0.2d, 0d, 0.4d, 1.7d, 3d, 6d, 12d, 25d, 60d })
                {
                    double x = mu + sigma * z;
                    double logPdf = ReferenceLogPDF(mu, sigma, g, x);
                    double logCdf = ReferenceLogTail(mu, sigma, g, x, false), logCcdf = ReferenceLogTail(mu, sigma, g, x, true);
                    Assert.AreEqual(logPdf, d.LogPDF(x), 0d);
                    Assert.AreEqual(Math.Exp(logPdf), d.PDF(x), 0d);
                    Assert.AreEqual(logCdf, d.LogCDF(x), 0d);
                    Assert.AreEqual(Math.Exp(logCdf), d.CDF(x), 0d);
                    Assert.AreEqual(logCcdf, d.LogCCDF(x), 0d);
                    Assert.AreEqual(Math.Exp(logCcdf), d.CCDF(x), 0d);
                }
            }
        }

        /// <summary>The gamma-path log density for skews outside the local-expansion region.</summary>
        /// <param name="mu">The mean.</param>
        /// <param name="sigma">The standard deviation.</param>
        /// <param name="g">The skew, zero or with absolute value at least 0.1.</param>
        /// <param name="x">The observation.</param>
        /// <returns>The log density.</returns>
        private static double ReferenceLogPDF(double mu, double sigma, double g, double x)
        {
            if (x < (g > 0d ? mu - sigma * (2d / g) : double.NegativeInfinity)
                || x > (g < 0d ? mu - sigma * (2d / g) : double.PositiveInfinity)
                || double.IsInfinity(x)) return double.NegativeInfinity;
            double z = DistributionNumerics.Standardize(x, mu, sigma);
            if (g == 0d) return -0.5d * z * z - Math.Log(sigma) - Math.Log(Tools.Sqrt2PI);
            double unit = ReferenceUnit(mu, sigma, g, x, z);
            return DistributionNumerics.GammaLogDensity(Math.Pow(2d / g, 2d), unit) - Math.Log(sigma) - Math.Log(Math.Abs(g) / 2d);
        }

        /// <summary>The gamma-path log tail for skews outside the local-expansion region.</summary>
        /// <param name="mu">The mean.</param>
        /// <param name="sigma">The standard deviation.</param>
        /// <param name="g">The skew, zero or with absolute value at least 0.1.</param>
        /// <param name="x">The observation.</param>
        /// <param name="upper"><see langword="true"/> for the survival probability.</param>
        /// <returns>The log tail probability.</returns>
        private static double ReferenceLogTail(double mu, double sigma, double g, double x, bool upper)
        {
            if (x <= (g > 0d ? mu - sigma * (2d / g) : double.NegativeInfinity)) return upper ? 0d : double.NegativeInfinity;
            if (x >= (g < 0d ? mu - sigma * (2d / g) : double.PositiveInfinity)) return upper ? double.NegativeInfinity : 0d;
            double z = DistributionNumerics.Standardize(x, mu, sigma);
            if (g == 0d) return upper ? DistributionNumerics.NormalLogSurvival(z) : DistributionNumerics.NormalLogCDF(z);
            double unit = ReferenceUnit(mu, sigma, g, x, z);
            double alpha = Math.Pow(2d / g, 2d);
            return upper == (g > 0d) ? DistributionNumerics.GammaLogSurvival(alpha, unit) : DistributionNumerics.GammaLogCDF(alpha, unit);
        }

        /// <summary>The unit gamma coordinate in the scaled-shape form.</summary>
        /// <param name="mu">The mean.</param>
        /// <param name="sigma">The standard deviation.</param>
        /// <param name="g">The nonzero skew.</param>
        /// <param name="x">The observation.</param>
        /// <param name="z">The standardized observation.</param>
        /// <returns>The unit gamma coordinate.</returns>
        private static double ReferenceUnit(double mu, double sigma, double g, double x, double z)
        {
            if (x == mu - sigma * (2d / g)) return 0d;
            return Math.Abs(g) < 1d ? Math.Pow(2d / g, 2d) * (1d + (g / 2d) * z)
                : DistributionNumerics.Standardize(x, mu - sigma * (2d / g), 0.5d * sigma * g);
        }
    }
}
