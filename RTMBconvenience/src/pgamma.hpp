/*
 *  Mathlib : A C Library of Special Functions
 *  Copyright (C) 1998 Ross Ihaka
 *  Copyright (C) 2000-2014 The R Core Team
 *  Copyright (C) 2003	    The R Foundation
 *
 *  This program is free software; you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation; either version 2 of the License, or
 *  (at your option) any later version.
 *
 *  This program is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with this program; if not, a copy is available at
 *  https://www.R-project.org/Licenses/
 *
 *  SYNOPSIS
 *
 *	double dnorm4(double x, double mu, double sigma, int give_log)
 *	      {dnorm (..) is synonymous and preferred inside R}
 *
 *  DESCRIPTION
 *
 *	Compute the density of the normal distribution.
 */

template<class Float>
Float dnorm4(Float x, Float mu, Float sigma, int give_log)
{
#ifdef IEEE_754
    if (ISNAN(x) || ISNAN(mu) || ISNAN(sigma))
	return x + mu + sigma;
#endif
    if (sigma < 0) ML_WARN_return_NAN;
    if(!R_FINITE(sigma)) return R_D__0;
    if(!R_FINITE(x) && mu == x) return ML_NAN;/* x-mu is NaN */
    if (sigma == 0) 
	return (x == mu) ? ML_POSINF : R_D__0;
    x = (x - mu) / sigma;

    if(!R_FINITE(x)) return R_D__0;

    x = fabs (x);
    if (x >= 2 * sqrt(DBL_MAX)) return R_D__0;
    if (give_log)
	return -(M_LN_SQRT_2PI + 0.5 * x * x + log(sigma));
    //  M_1_SQRT_2PI = 1 / sqrt(2 * pi)
#ifdef MATHLIB_FAST_dnorm
    // and for R <= 3.0.x and R-devel upto 2014-01-01:
    return M_1_SQRT_2PI * exp(-0.5 * x * x) / sigma;
#else
    // more accurate, less fast :
    if (x < 5)    return M_1_SQRT_2PI * exp(-0.5 * x * x) / sigma;

    /* ELSE:

     * x*x  may lose upto about two digits accuracy for "large" x
     * Morten Welinder's proposal for PR#15620
     * https://bugs.r-project.org/show_bug.cgi?id=15620

     * -- 1 --  No hoop jumping when we underflow to zero anyway:

     *  -x^2/2 <         log(2)*.Machine$double.min.exp  <==>
     *     x   > sqrt(-2*log(2)*.Machine$double.min.exp) =IEEE= 37.64031
     * but "thanks" to denormalized numbers, underflow happens a bit later,
     *  effective.D.MIN.EXP <- with(.Machine, double.min.exp + double.ulp.digits)
     * for IEEE, DBL_MIN_EXP is -1022 but "effective" is -1074
     * ==> boundary = sqrt(-2*log(2)*(.Machine$double.min.exp + .Machine$double.ulp.digits))
     *              =IEEE=  38.58601
     * [on one x86_64 platform, effective boundary a bit lower: 38.56804]
     */
    if (x > sqrt(-2*M_LN2*(DBL_MIN_EXP + 1-DBL_MANT_DIG))) return 0.;

    /* Now, to get full accuracy, split x into two parts,
     *  x = x1+x2, such that |x2| <= 2^-16.
     * Assuming that we are using IEEE doubles, that means that
     * x1*x1 is error free for x<1024 (but we have x < 38.6 anyway).

     * If we do not have IEEE this is still an improvement over the naive formula.
     */
    Float x1 = //  R_forceint(x * 65536) / 65536 =
      ldexp2( (Float) (int)trunc((ldexp2((Float)x, (Float)16))+0.5), (Float) -16);
    Float x2 = x - x1;
    return M_1_SQRT_2PI / sigma *
	(exp(-0.5 * x1 * x1) * exp( (-0.5 * x2 - x1) * x2 ) );
#endif
}


template<class Float>
Float log1pmx (Float x);

 template<class Float>
 Float lgamma1p (Float a);

static const float bd0_scale[128 + 1][4] = {
	{ +0x1.62e430p-1, -0x1.05c610p-29, -0x1.950d88p-54, +0x1.d9cc02p-79 }, /* 128: log(2048/1024.) */
	{ +0x1.5ee02cp-1, -0x1.6dbe98p-25, -0x1.51e540p-50, +0x1.2bfa48p-74 }, /* 129: log(2032/1024.) */
	{ +0x1.5ad404p-1, +0x1.86b3e4p-26, +0x1.9f6534p-50, +0x1.54be04p-74 }, /* 130: log(2016/1024.) */
	{ +0x1.570124p-1, -0x1.9ed750p-25, -0x1.f37dd0p-51, +0x1.10b770p-77 }, /* 131: log(2001/1024.) */
	{ +0x1.5326e4p-1, -0x1.9b9874p-25, -0x1.378194p-49, +0x1.56feb2p-74 }, /* 132: log(1986/1024.) */
	{ +0x1.4f4528p-1, +0x1.aca70cp-28, +0x1.103e74p-53, +0x1.9c410ap-81 }, /* 133: log(1971/1024.) */
	{ +0x1.4b5bd8p-1, -0x1.6a91d8p-25, -0x1.8e43d0p-50, -0x1.afba9ep-77 }, /* 134: log(1956/1024.) */
	{ +0x1.47ae54p-1, -0x1.abb51cp-25, +0x1.19b798p-51, +0x1.45e09cp-76 }, /* 135: log(1942/1024.) */
	{ +0x1.43fa00p-1, -0x1.d06318p-25, -0x1.8858d8p-49, -0x1.1927c4p-75 }, /* 136: log(1928/1024.) */
	{ +0x1.3ffa40p-1, +0x1.1a427cp-25, +0x1.151640p-53, -0x1.4f5606p-77 }, /* 137: log(1913/1024.) */
	{ +0x1.3c7c80p-1, -0x1.19bf48p-34, +0x1.05fc94p-58, -0x1.c096fcp-82 }, /* 138: log(1900/1024.) */
	{ +0x1.38b320p-1, +0x1.6b5778p-25, +0x1.be38d0p-50, -0x1.075e96p-74 }, /* 139: log(1886/1024.) */
	{ +0x1.34e288p-1, +0x1.d9ce1cp-25, +0x1.316eb8p-49, +0x1.2d885cp-73 }, /* 140: log(1872/1024.) */
	{ +0x1.315124p-1, +0x1.c2fc60p-29, -0x1.4396fcp-53, +0x1.acf376p-78 }, /* 141: log(1859/1024.) */
	{ +0x1.2db954p-1, +0x1.720de4p-25, -0x1.d39b04p-49, -0x1.f11176p-76 }, /* 142: log(1846/1024.) */
	{ +0x1.2a1b08p-1, -0x1.562494p-25, +0x1.a7863cp-49, +0x1.85dd64p-73 }, /* 143: log(1833/1024.) */
	{ +0x1.267620p-1, +0x1.3430e0p-29, -0x1.96a958p-56, +0x1.f8e636p-82 }, /* 144: log(1820/1024.) */
	{ +0x1.23130cp-1, +0x1.7bebf4p-25, +0x1.416f1cp-52, -0x1.78dd36p-77 }, /* 145: log(1808/1024.) */
	{ +0x1.1faa34p-1, +0x1.70e128p-26, +0x1.81817cp-50, -0x1.c2179cp-76 }, /* 146: log(1796/1024.) */
	{ +0x1.1bf204p-1, +0x1.3a9620p-28, +0x1.2f94c0p-52, +0x1.9096c0p-76 }, /* 147: log(1783/1024.) */
	{ +0x1.187ce4p-1, -0x1.077870p-27, +0x1.655a80p-51, +0x1.eaafd6p-78 }, /* 148: log(1771/1024.) */
	{ +0x1.1501c0p-1, -0x1.406cacp-25, -0x1.e72290p-49, +0x1.5dd800p-73 }, /* 149: log(1759/1024.) */
	{ +0x1.11cb80p-1, +0x1.787cd0p-25, -0x1.efdc78p-51, -0x1.5380cep-77 }, /* 150: log(1748/1024.) */
	{ +0x1.0e4498p-1, +0x1.747324p-27, -0x1.024548p-51, +0x1.77a5a6p-75 }, /* 151: log(1736/1024.) */
	{ +0x1.0b036cp-1, +0x1.690c74p-25, +0x1.5d0cc4p-50, -0x1.c0e23cp-76 }, /* 152: log(1725/1024.) */
	{ +0x1.077070p-1, -0x1.a769bcp-27, +0x1.452234p-52, +0x1.6ba668p-76 }, /* 153: log(1713/1024.) */
	{ +0x1.04240cp-1, -0x1.a686acp-27, -0x1.ef46b0p-52, -0x1.5ce10cp-76 }, /* 154: log(1702/1024.) */
	{ +0x1.00d22cp-1, +0x1.fc0e10p-25, +0x1.6ee034p-50, -0x1.19a2ccp-74 }, /* 155: log(1691/1024.) */
	{ +0x1.faf588p-2, +0x1.ef1e64p-27, -0x1.26504cp-54, -0x1.b15792p-82 }, /* 156: log(1680/1024.) */
	{ +0x1.f4d87cp-2, +0x1.d7b980p-26, -0x1.a114d8p-50, +0x1.9758c6p-75 }, /* 157: log(1670/1024.) */
	{ +0x1.ee1414p-2, +0x1.2ec060p-26, +0x1.dc00fcp-52, +0x1.f8833cp-76 }, /* 158: log(1659/1024.) */
	{ +0x1.e7e32cp-2, -0x1.ac796cp-27, -0x1.a68818p-54, +0x1.235d02p-78 }, /* 159: log(1649/1024.) */
	{ +0x1.e108a0p-2, -0x1.768ba4p-28, -0x1.f050a8p-52, +0x1.00d632p-82 }, /* 160: log(1638/1024.) */
	{ +0x1.dac354p-2, -0x1.d3a6acp-30, +0x1.18734cp-57, -0x1.f97902p-83 }, /* 161: log(1628/1024.) */
	{ +0x1.d47424p-2, +0x1.7dbbacp-31, -0x1.d5ada4p-56, +0x1.56fcaap-81 }, /* 162: log(1618/1024.) */
	{ +0x1.ce1af0p-2, +0x1.70be7cp-27, +0x1.6f6fa4p-51, +0x1.7955a2p-75 }, /* 163: log(1608/1024.) */
	{ +0x1.c7b798p-2, +0x1.ec36ecp-26, -0x1.07e294p-50, -0x1.ca183cp-75 }, /* 164: log(1598/1024.) */
	{ +0x1.c1ef04p-2, +0x1.c1dfd4p-26, +0x1.888eecp-50, -0x1.fd6b86p-75 }, /* 165: log(1589/1024.) */
	{ +0x1.bb7810p-2, +0x1.478bfcp-26, +0x1.245b8cp-50, +0x1.ea9d52p-74 }, /* 166: log(1579/1024.) */
	{ +0x1.b59da0p-2, -0x1.882b08p-27, +0x1.31573cp-53, -0x1.8c249ap-77 }, /* 167: log(1570/1024.) */
	{ +0x1.af1294p-2, -0x1.b710f4p-27, +0x1.622670p-51, +0x1.128578p-76 }, /* 168: log(1560/1024.) */
	{ +0x1.a925d4p-2, -0x1.0ae750p-27, +0x1.574ed4p-51, +0x1.084996p-75 }, /* 169: log(1551/1024.) */
	{ +0x1.a33040p-2, +0x1.027d30p-29, +0x1.b9a550p-53, -0x1.b2e38ap-78 }, /* 170: log(1542/1024.) */
	{ +0x1.9d31c0p-2, -0x1.5ec12cp-26, -0x1.5245e0p-52, +0x1.2522d0p-79 }, /* 171: log(1533/1024.) */
	{ +0x1.972a34p-2, +0x1.135158p-30, +0x1.a5c09cp-56, +0x1.24b70ep-80 }, /* 172: log(1524/1024.) */
	{ +0x1.911984p-2, +0x1.0995d4p-26, +0x1.3bfb5cp-50, +0x1.2c9dd6p-75 }, /* 173: log(1515/1024.) */
	{ +0x1.8bad98p-2, -0x1.1d6144p-29, +0x1.5b9208p-53, +0x1.1ec158p-77 }, /* 174: log(1507/1024.) */
	{ +0x1.858b58p-2, -0x1.1b4678p-27, +0x1.56cab4p-53, -0x1.2fdc0cp-78 }, /* 175: log(1498/1024.) */
	{ +0x1.7f5fa0p-2, +0x1.3aaf48p-27, +0x1.461964p-51, +0x1.4ae476p-75 }, /* 176: log(1489/1024.) */
	{ +0x1.79db68p-2, -0x1.7e5054p-26, +0x1.673750p-51, -0x1.a11f7ap-76 }, /* 177: log(1481/1024.) */
	{ +0x1.744f88p-2, -0x1.cc0e18p-26, -0x1.1e9d18p-50, -0x1.6c06bcp-78 }, /* 178: log(1473/1024.) */
	{ +0x1.6e08ecp-2, -0x1.5d45e0p-26, -0x1.c73ec8p-50, +0x1.318d72p-74 }, /* 179: log(1464/1024.) */
	{ +0x1.686c80p-2, +0x1.e9b14cp-26, -0x1.13bbd4p-50, -0x1.efeb1cp-78 }, /* 180: log(1456/1024.) */
	{ +0x1.62c830p-2, -0x1.a8c70cp-27, -0x1.5a1214p-51, -0x1.bab3fcp-79 }, /* 181: log(1448/1024.) */
	{ +0x1.5d1bdcp-2, -0x1.4fec6cp-31, +0x1.423638p-56, +0x1.ee3feep-83 }, /* 182: log(1440/1024.) */
	{ +0x1.576770p-2, +0x1.7455a8p-26, -0x1.3ab654p-50, -0x1.26be4cp-75 }, /* 183: log(1432/1024.) */
	{ +0x1.5262e0p-2, -0x1.146778p-26, -0x1.b9f708p-52, -0x1.294018p-77 }, /* 184: log(1425/1024.) */
	{ +0x1.4c9f08p-2, +0x1.e152c4p-26, -0x1.dde710p-53, +0x1.fd2208p-77 }, /* 185: log(1417/1024.) */
	{ +0x1.46d2d8p-2, +0x1.c28058p-26, -0x1.936284p-50, +0x1.9fdd68p-74 }, /* 186: log(1409/1024.) */
	{ +0x1.41b940p-2, +0x1.cce0c0p-26, -0x1.1a4050p-50, +0x1.bc0376p-76 }, /* 187: log(1402/1024.) */
	{ +0x1.3bdd24p-2, +0x1.d6296cp-27, +0x1.425b48p-51, -0x1.cddb2cp-77 }, /* 188: log(1394/1024.) */
	{ +0x1.36b578p-2, -0x1.287ddcp-27, -0x1.2d0f4cp-51, +0x1.38447ep-75 }, /* 189: log(1387/1024.) */
	{ +0x1.31871cp-2, +0x1.2a8830p-27, +0x1.3eae54p-52, -0x1.898136p-77 }, /* 190: log(1380/1024.) */
	{ +0x1.2b9304p-2, -0x1.51d8b8p-28, +0x1.27694cp-52, -0x1.fd852ap-76 }, /* 191: log(1372/1024.) */
	{ +0x1.265620p-2, -0x1.d98f3cp-27, +0x1.a44338p-51, -0x1.56e85ep-78 }, /* 192: log(1365/1024.) */
	{ +0x1.211254p-2, +0x1.986160p-26, +0x1.73c5d0p-51, +0x1.4a861ep-75 }, /* 193: log(1358/1024.) */
	{ +0x1.1bc794p-2, +0x1.fa3918p-27, +0x1.879c5cp-51, +0x1.16107cp-78 }, /* 194: log(1351/1024.) */
	{ +0x1.1675ccp-2, -0x1.4545a0p-26, +0x1.c07398p-51, +0x1.f55c42p-76 }, /* 195: log(1344/1024.) */
	{ +0x1.111ce4p-2, +0x1.f72670p-37, -0x1.b84b5cp-61, +0x1.a4a4dcp-85 }, /* 196: log(1337/1024.) */
	{ +0x1.0c81d4p-2, +0x1.0c150cp-27, +0x1.218600p-51, -0x1.d17312p-76 }, /* 197: log(1331/1024.) */
	{ +0x1.071b84p-2, +0x1.fcd590p-26, +0x1.a3a2e0p-51, +0x1.fe5ef8p-76 }, /* 198: log(1324/1024.) */
	{ +0x1.01ade4p-2, -0x1.bb1844p-28, +0x1.db3cccp-52, +0x1.1f56fcp-77 }, /* 199: log(1317/1024.) */
	{ +0x1.fa01c4p-3, -0x1.12a0d0p-29, -0x1.f71fb0p-54, +0x1.e287a4p-78 }, /* 200: log(1311/1024.) */
	{ +0x1.ef0adcp-3, +0x1.7b8b28p-28, -0x1.35bce4p-52, -0x1.abc8f8p-79 }, /* 201: log(1304/1024.) */
	{ +0x1.e598ecp-3, +0x1.5a87e4p-27, -0x1.134bd0p-51, +0x1.c2cebep-76 }, /* 202: log(1298/1024.) */
	{ +0x1.da85d8p-3, -0x1.df31b0p-27, +0x1.94c16cp-57, +0x1.8fd7eap-82 }, /* 203: log(1291/1024.) */
	{ +0x1.d0fb80p-3, -0x1.bb5434p-28, -0x1.ea5640p-52, -0x1.8ceca4p-77 }, /* 204: log(1285/1024.) */
	{ +0x1.c765b8p-3, +0x1.e4d68cp-27, +0x1.5b59b4p-51, +0x1.76f6c4p-76 }, /* 205: log(1279/1024.) */
	{ +0x1.bdc46cp-3, -0x1.1cbb50p-27, +0x1.2da010p-51, +0x1.eb282cp-75 }, /* 206: log(1273/1024.) */
	{ +0x1.b27980p-3, -0x1.1b9ce0p-27, +0x1.7756f8p-52, +0x1.2ff572p-76 }, /* 207: log(1266/1024.) */
	{ +0x1.a8bed0p-3, -0x1.bbe874p-30, +0x1.85cf20p-56, +0x1.b9cf18p-80 }, /* 208: log(1260/1024.) */
	{ +0x1.9ef83cp-3, +0x1.2769a4p-27, -0x1.85bda0p-52, +0x1.8c8018p-79 }, /* 209: log(1254/1024.) */
	{ +0x1.9525a8p-3, +0x1.cf456cp-27, -0x1.7137d8p-52, -0x1.f158e8p-76 }, /* 210: log(1248/1024.) */
	{ +0x1.8b46f8p-3, +0x1.11b12cp-30, +0x1.9f2104p-54, -0x1.22836ep-78 }, /* 211: log(1242/1024.) */
	{ +0x1.83040cp-3, +0x1.2379e4p-28, +0x1.b71c70p-52, -0x1.990cdep-76 }, /* 212: log(1237/1024.) */
	{ +0x1.790ed4p-3, +0x1.dc4c68p-28, -0x1.910ac8p-52, +0x1.dd1bd6p-76 }, /* 213: log(1231/1024.) */
	{ +0x1.6f0d28p-3, +0x1.5cad68p-28, +0x1.737c94p-52, -0x1.9184bap-77 }, /* 214: log(1225/1024.) */
	{ +0x1.64fee8p-3, +0x1.04bf88p-28, +0x1.6fca28p-52, +0x1.8884a8p-76 }, /* 215: log(1219/1024.) */
	{ +0x1.5c9400p-3, +0x1.d65cb0p-29, -0x1.b2919cp-53, +0x1.b99bcep-77 }, /* 216: log(1214/1024.) */
	{ +0x1.526e60p-3, -0x1.c5e4bcp-27, -0x1.0ba380p-52, +0x1.d6e3ccp-79 }, /* 217: log(1208/1024.) */
	{ +0x1.483bccp-3, +0x1.9cdc7cp-28, -0x1.5ad8dcp-54, -0x1.392d3cp-83 }, /* 218: log(1202/1024.) */
	{ +0x1.3fb25cp-3, -0x1.a6ad74p-27, +0x1.5be6b4p-52, -0x1.4e0114p-77 }, /* 219: log(1197/1024.) */
	{ +0x1.371fc4p-3, -0x1.fe1708p-27, -0x1.78864cp-52, -0x1.27543ap-76 }, /* 220: log(1192/1024.) */
	{ +0x1.2cca10p-3, -0x1.4141b4p-28, -0x1.ef191cp-52, +0x1.00ee08p-76 }, /* 221: log(1186/1024.) */
	{ +0x1.242310p-3, +0x1.3ba510p-27, -0x1.d003c8p-51, +0x1.162640p-76 }, /* 222: log(1181/1024.) */
	{ +0x1.1b72acp-3, +0x1.52f67cp-27, -0x1.fd6fa0p-51, +0x1.1a3966p-77 }, /* 223: log(1176/1024.) */
	{ +0x1.10f8e4p-3, +0x1.129cd8p-30, +0x1.31ef30p-55, +0x1.a73e38p-79 }, /* 224: log(1170/1024.) */
	{ +0x1.08338cp-3, -0x1.005d7cp-27, -0x1.661a9cp-51, +0x1.1f138ap-79 }, /* 225: log(1165/1024.) */
	{ +0x1.fec914p-4, -0x1.c482a8p-29, -0x1.55746cp-54, +0x1.99f932p-80 }, /* 226: log(1160/1024.) */
	{ +0x1.ed1794p-4, +0x1.d06f00p-29, +0x1.75e45cp-53, -0x1.d0483ep-78 }, /* 227: log(1155/1024.) */
	{ +0x1.db5270p-4, +0x1.87d928p-32, -0x1.0f52a4p-57, +0x1.81f4a6p-84 }, /* 228: log(1150/1024.) */
	{ +0x1.c97978p-4, +0x1.af1d24p-29, -0x1.0977d0p-60, -0x1.8839d0p-84 }, /* 229: log(1145/1024.) */
	{ +0x1.b78c84p-4, -0x1.44f124p-28, -0x1.ef7bc4p-52, +0x1.9e0650p-78 }, /* 230: log(1140/1024.) */
	{ +0x1.a58b60p-4, +0x1.856464p-29, +0x1.c651d0p-55, +0x1.b06b0cp-79 }, /* 231: log(1135/1024.) */
	{ +0x1.9375e4p-4, +0x1.5595ecp-28, +0x1.dc3738p-52, +0x1.86c89ap-81 }, /* 232: log(1130/1024.) */
	{ +0x1.814be4p-4, -0x1.c073fcp-28, -0x1.371f88p-53, -0x1.5f4080p-77 }, /* 233: log(1125/1024.) */
	{ +0x1.6f0d28p-4, +0x1.5cad68p-29, +0x1.737c94p-53, -0x1.9184bap-78 }, /* 234: log(1120/1024.) */
	{ +0x1.60658cp-4, -0x1.6c8af4p-28, +0x1.d8ef74p-55, +0x1.c4f792p-80 }, /* 235: log(1116/1024.) */
	{ +0x1.4e0110p-4, +0x1.146b5cp-29, +0x1.73f7ccp-54, -0x1.d28db8p-79 }, /* 236: log(1111/1024.) */
	{ +0x1.3b8758p-4, +0x1.8b1b70p-28, -0x1.20aca4p-52, -0x1.651894p-76 }, /* 237: log(1106/1024.) */
	{ +0x1.28f834p-4, +0x1.43b6a4p-30, -0x1.452af8p-55, +0x1.976892p-80 }, /* 238: log(1101/1024.) */
	{ +0x1.1a0fbcp-4, -0x1.e4075cp-28, +0x1.1fe618p-52, +0x1.9d6dc2p-77 }, /* 239: log(1097/1024.) */
	{ +0x1.075984p-4, -0x1.4ce370p-29, -0x1.d9fc98p-53, +0x1.4ccf12p-77 }, /* 240: log(1092/1024.) */
	{ +0x1.f0a30cp-5, +0x1.162a68p-37, -0x1.e83368p-61, -0x1.d222a6p-86 }, /* 241: log(1088/1024.) */
	{ +0x1.cae730p-5, -0x1.1a8f7cp-31, -0x1.5f9014p-55, +0x1.2720c0p-79 }, /* 242: log(1083/1024.) */
	{ +0x1.ac9724p-5, -0x1.e8ee08p-29, +0x1.a7de04p-54, -0x1.9bba74p-78 }, /* 243: log(1079/1024.) */
	{ +0x1.868a84p-5, -0x1.ef8128p-30, +0x1.dc5eccp-54, -0x1.58d250p-79 }, /* 244: log(1074/1024.) */
	{ +0x1.67f950p-5, -0x1.ed684cp-30, -0x1.f060c0p-55, -0x1.b1294cp-80 }, /* 245: log(1070/1024.) */
	{ +0x1.494accp-5, +0x1.a6c890p-32, -0x1.c3ad48p-56, -0x1.6dc66cp-84 }, /* 246: log(1066/1024.) */
	{ +0x1.22c71cp-5, -0x1.8abe2cp-32, -0x1.7e7078p-56, -0x1.ddc3dcp-86 }, /* 247: log(1061/1024.) */
	{ +0x1.03d5d8p-5, +0x1.79cfbcp-31, -0x1.da7c4cp-58, +0x1.4e7582p-83 }, /* 248: log(1057/1024.) */
	{ +0x1.c98d18p-6, +0x1.a01904p-31, -0x1.854164p-55, +0x1.883c36p-79 }, /* 249: log(1053/1024.) */
	{ +0x1.8b31fcp-6, -0x1.356500p-30, +0x1.c3ab48p-55, +0x1.b69bdap-80 }, /* 250: log(1049/1024.) */
	{ +0x1.3cea44p-6, +0x1.a352bcp-33, -0x1.8865acp-57, -0x1.48159cp-81 }, /* 251: log(1044/1024.) */
	{ +0x1.fc0a8cp-7, -0x1.e07f84p-32, +0x1.e7cf6cp-58, +0x1.3a69c0p-82 }, /* 252: log(1040/1024.) */
	{ +0x1.7dc474p-7, +0x1.f810a8p-31, -0x1.245b5cp-56, -0x1.a1f4f8p-80 }, /* 253: log(1036/1024.) */
	{ +0x1.fe02a8p-8, -0x1.4ef988p-32, +0x1.1f86ecp-57, +0x1.20723cp-81 }, /* 254: log(1032/1024.) */
	{ +0x1.ff00acp-9, -0x1.d4ef44p-33, +0x1.2821acp-63, +0x1.5a6d32p-87 }, /* 255: log(1028/1024.) */
	{ 0, 0, 0, 0 } /* log(1024/1024) = log(1) = 0 */
};


/*
 * Compute x * log (x / M) + (M - x)
 * aka -x * log1pmx ((M - x) / x)
 *
 * Deliver the result back in two parts, *yh and *yl.
 */
template<class Float>
void attribute_hidden ebd0(Float x, Float M, Float *yh, Float *yl)
{
	const int Sb = 10;
	const Float S = 1u << Sb; // = 2^10 = 1024
	const int N = 128; // == ? == G_N_ELEMENTS(bd0_scale) - 1; <<<< FIXME:

	*yl = *yh = 0;

	if (x == M) return;
	if (x == 0) { *yh = M;         return; }
	if (M == 0) { *yh = ML_POSINF; return; }

	if (M/x == ML_POSINF) { *yh = M; return; }//  as when (x == 0)

	int e;
	// NB: M/x overflow handled above; underflow should be handled by fg = Inf
	Float r = frexp2 ((Float)(M / x), &e); // => r in  [0.5, 1) and 'e' (int) such that  M/x = r * 2^e

	// prevent later overflow
	if (M_LN2 * (-e)  > 1. + DBL_MAX / x) { *yh = ML_POSINF; return; }

	int i = (int) floor ((r - 0.5) * (2 * N) + 0.5);
	// now,  0 <= i <= N
	Float f = trunc (S / (0.5 + i / (2.0 * N)) + 0.5);
	Float fg = ldexp2 ((Float)f, (Float)-(e + Sb)); // ldexp(f, E) := f * 2^E
#ifdef DEBUG_bd0
	REprintf("ebd0(x=%g, M=%g): M/x = (r=%.15g) * 2^(e=%d); i=%d,\n  f=%g, fg=f*2^-(e+%d)=%g\n",
		 x, M, r,e, i, f, Sb, fg);
	if (fg == ML_POSINF) {
	    REprintf(" --> fg = +Inf --> return( +Inf )\n");
	    *yh = fg; return;
	}
	REprintf("     bd0_sc[0][0..3]= ("); for(int j=0; j < 4; j++) REprintf("%g ", bd0_scale[0][j]); REprintf(")\n");
	REprintf("i -> bd0_sc[i][0..3]= ("); for(int j=0; j < 4; j++) REprintf("%g ", bd0_scale[i][j]); REprintf(")\n");
	REprintf( "  small(?)  (M*fg-x)/x = (M*fg)/x - 1 = %.16g\n", (M*fg-x)/x);
#else
	if (fg == ML_POSINF) {
	    *yh = fg; return;
	}
#endif
	/* We now have (M * fg / x) close to 1.  */

	/*
	 * We need to compute this:
	 * (x/M)^x * exp(M-x) =
	 * (M/x)^-x * exp(M-x) =
	 * (M*fg/x)^-x * (fg)^x * exp(M-x) =
	 * (M*fg/x)^-x * (fg)^x * exp(M*fg-x) * exp(M-M*fg)
	 *
	 * In log terms:
	 * log((x/M)^x * exp(M-x)) =
	 * log((M*fg/x)^-x * (fg)^x * exp(M*fg-x) * exp(M-M*fg)) =
	 * log((M*fg/x)^-x * exp(M*fg-x)) + x*log(fg) + (M-M*fg) =
	 * -x*log1pmx((M*fg-x)/x) + x*log(fg) + M - M*fg =
	 *
	 * Note, that fg has at most 10 bits.  If M and x are suitably
	 * "nice" -- such as being integers or half-integers -- then
	 * we can compute M*fg as well as x * bd0_scale[.][.] without
	 * rounding errors.
	 */

#define ADD1(d_) do {				\
   Float d = (d_);			\
   Float dx = d + 0.5;			\
   Float d1 = floor ((Float)dx);	\
	    Float d2 = d - d1;/* in [-.5,.5) */ \
	    *yh += d1;				\
	    *yl += d2;				\
	} while(0)

#ifdef DEBUG_bd0
	{
	    Float log1__ = log1pmx((M * fg - x) / x),
		xl = -x * log1__;
	    REprintf(" 1a. before adding  -x * log1pmx(.) = -x * %g = %g\n", log1__, xl);
	    ADD1(xl);
	    REprintf(" 1. after A.(-x*l..):       yl,yh = (%13g, %13g); yl+yh= %g\n",
		     *yl, *yh, (*yl)+(*yh));
	}
        if(fg == 1) {
            REprintf("___ fg = 1 ___ skipping further steps\n");
            return;
        }
	// else  [ fg != 1 ]
	REprintf(" 2:  A(x*b[i,j]) and A(-x*e*b[0,j]), j=1:4:\n");
	for (int j = 0; j < 4; j++) {
 	    ADD1( x * bd0_scale[i][j]);     // handles  x*log(fg*2^e)
	    REprintf(" j=%d: (%13g, %13g);", j, *yl, *yh);
	    ADD1(-x * bd0_scale[0][j] * e); // handles  x*log(1/ 2^e)
	    REprintf(" (%13g, %13g); yl+yh= %g\n", *yl, *yh, (*yl)+(*yh));
            if(!R_FINITE(*yh)) {
                REprintf(" non-finite yh --> return((yh=Inf, yl=0))\n");
		*yh = ML_POSINF; *yl = 0; return;
            }
	}
#else
	ADD1(-x * log1pmx ((M * fg - x) / x));
        if(fg == 1) return;
	// else (fg != 1) :
	for (int j = 0; j < 4; j++) {
	    ADD1( x * bd0_scale[i][j]);     // handles  x*log(fg*2^e)
	    ADD1(-x * bd0_scale[0][j] * e); // handles  x*log(1/ 2^e)
	    //                        ^^^ at end prevents overflow in  ebd0(1e307, 1e300)
            if(!R_FINITE(*yh)) { *yh = ML_POSINF; *yl = 0; return; }
	}
#endif

	ADD1(M);
#ifdef DEBUG_bd0
	REprintf(" 3. after ADD1(M):            yl,yh = (%13g, %13g); yl+yh= %g\n", *yl, *yh, (*yl)+(*yh));
#endif
	ADD1(-M * fg);
#ifdef DEBUG_bd0
	REprintf(" 4. after ADD1(- M*fg):       yl,yh = (%13g, %13g); yl+yh= %g\n\n", *yl, *yh, (*yl)+(*yh));
#endif
}

#undef ADD1




#define M_SQRT_2PI	2.50662827463100050241576528481104525301  /* sqrt(2*pi) */
#define x_LRG           2.86111748575702815380240589208115399625e+307 /* = 2^1023 / pi */

template<class Float>
Float _RORIGINAL_dpois_raw(Float x, Float lambda, int give_log)
{
    /*       x >= 0 ; integer for dpois(), but not e.g. for pgamma()!
        lambda >= 0
    */
    if (lambda == 0) return( (x == 0) ? R_D__1 : R_D__0 );
    if (!R_FINITE(lambda)) return R_D__0; // including for the case where  x = lambda = +Inf
    if (x < 0) return( R_D__0 );
    if (x <= lambda * DBL_MIN) return(R_D_exp(-lambda) );
    if (lambda < x * DBL_MIN) {
	if (!R_FINITE(x)) // lambda < x = +Inf
	    return R_D__0;
	// else
	return(R_D_exp(-lambda + x*log(lambda) -lgammafn(x+1)));
    }
    // R <= 4.0.x  had   return(R_D_fexp( M_2PI*x, -stirlerr(x)-bd0(x,lambda) ));
    Float yh, yl;
    ebd0 (x, lambda, &yh, &yl);
    yl += stirlerr(x);
    bool Lrg_x = (x >= x_LRG); //really large x  <==>  2*pi*x  overflows
    Float r = Lrg_x
	? M_SQRT_2PI * sqrt(x) // sqrt(.): avoid overflow for very large x
	: M_2PI * x;
    return give_log
	? -yl - yh - (Lrg_x ? log(r) : 0.5 * log(r))
	: exp(-yl) * exp(-yh) / (Lrg_x ? r : sqrt(r));
}

// The original R version above gave 0 when it was not supposed to - presumably from a change to make it work with AD/Float types
template<class Float>
Float dpois_raw(Float x, Float lambda, int give_log){
  if (lambda == 0) return( (x == 0) ? R_D__1 : R_D__0 );
  if (!R_FINITE(lambda)) return R_D__0; // including for the case where  x = lambda = +Inf
  if (x < 0) return( R_D__0 );
  if (x <= lambda * DBL_MIN) return(R_D_exp(-lambda) );
  if (lambda < x * DBL_MIN) {
    if (!R_FINITE(x)) // lambda < x = +Inf
      return R_D__0;
  }
  Float v= -lambda + x*log(lambda) -lgammafn(x+1);
  if(give_log)
    return v;
  return exp(v);
}


/*
 *  Mathlib : A C Library of Special Functions
 *  Copyright (C) 2006-2019 The R Core Team
 *  Copyright (C) 2005-6 Morten Welinder <terra@gnome.org>
 *  Copyright (C) 2005-10 The R Foundation
 *
 *  This program is free software; you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation; either version 2 of the License, or
 *  (at your option) any later version.
 *
 *  This program is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with this program; if not, a copy is available at
 *  https://www.R-project.org/Licenses/
 *
 *  SYNOPSIS
 *
 *	#include <Rmath.h>
 *
 *	double pgamma (double x, double alph, double scale,
 *		       int lower_tail, int log_p)
 *
 *	double log1pmx	(double x)
 *	double lgamma1p (double a)
 *
 *	double logspace_add (double logx, double logy)
 *	double logspace_sub (double logx, double logy)
 *	double logspace_sum (double* logx, int n)
 *
 *
 *  DESCRIPTION
 *
 *	This function computes the distribution function for the
 *	gamma distribution with shape parameter alph and scale parameter
 *	scale.	This is also known as the incomplete gamma function.
 *	See Abramowitz and Stegun (6.5.1) for example.
 *
 *  NOTES
 *
 *	Complete redesign by Morten Welinder, originally for Gnumeric.
 *	Improvements (e.g. "while NEEDED_SCALE") by Martin Maechler
 *
 *  REFERENCES
 *
 */

/*
 * Modified 2024-09-10 Christoffer Moesgaard Albertsen, Technical University of Denmark (cmoe@aqua.dtu.dk)
 */
    
/*----------- DEBUGGING -------------
 * make CFLAGS='-DDEBUG_p -g'
 * (cd `R-devel RHOME`/src/nmath; gcc -I. -I../../src/include -I../../../R/src/include  -DHAVE_CONFIG_H -fopenmp -DDEBUG_p -g -c ../../../R/src/nmath/pgamma.c -o pgamma.o)
 */

/* Scalefactor:= (2^32)^8 = 2^256 = 1.157921e+77 */
#define SQR(x) ((x)*(x))
static const double scalefactor = SQR(SQR(SQR(4294967296.0)));
#undef SQR

/* If |x| > |k| * M_cutoff,  then  log[ exp(-x) * k^x ]	 =~=  -x */
static const double M_cutoff = M_LN2 * DBL_MAX_EXP / DBL_EPSILON;/*=3.196577e18*/

/* Continued fraction for calculation of
 *    1/i + x/(i+d) + x^2/(i+2*d) + x^3/(i+3*d) + ... = sum_{k=0}^Inf x^k/(i+k*d)
 *
 * auxiliary in log1pmx() and lgamma1p()
 */
    template<class Float>
static Float
logcf (Float x, Float i, Float d,
       Float eps /* ~ relative tolerance */)
{
    Float c1 = 2 * d;
    Float c2 = i + d;
    Float c4 = c2 + d;
    Float a1 = c2;
    Float b1 = i * (c2 - i * x);
    Float b2 = d * d * x;
    Float a2 = c4 * c2 - b2;

#if 0
    assert (i > 0);
    assert (d >= 0);
#endif

    b2 = c4 * b1 - i * b2;

    while (fabs(a2 * b1 - a1 * b2) > fabs(eps * b1 * b2)) {
	Float c3 = c2*c2*x;
	c2 += d;
	c4 += d;
	a1 = c4 * a2 - c3 * a1;
	b1 = c4 * b2 - c3 * b1;

	c3 = c1 * c1 * x;
	c1 += d;
	c4 += d;
	a2 = c4 * a1 - c3 * a2;
	b2 = c4 * b1 - c3 * b2;

	if (fabs (b2) > scalefactor) {
	    a1 /= scalefactor;
	    b1 /= scalefactor;
	    a2 /= scalefactor;
	    b2 /= scalefactor;
	} else if (fabs (b2) < 1 / scalefactor) {
	    a1 *= scalefactor;
	    b1 *= scalefactor;
	    a2 *= scalefactor;
	    b2 *= scalefactor;
	}
    }

    return a2 / b2;
}

/* Accurate calculation of log(1+x)-x, particularly for small x.  */
    template<class Float>
    Float log1pmx (Float x)
{
    static const Float minLog1Value = -0.79149064;

    if (x > 1 || x < minLog1Value)
	return log1p(x) - x;
    else { /* -.791 <=  x <= 1  -- expand in  [x/(2+x)]^2 =: y :
	    * log(1+x) - x =  x/(2+x) * [ 2 * y * S(y) - x],  with
	    * ---------------------------------------------
	    * S(y) = 1/3 + y/5 + y^2/7 + ... = \sum_{k=0}^\infty  y^k / (2k + 3)
	   */
	Float r = x / (2 + x), y = r * r;
	if (fabs(x) < 1e-2) {
	    static const Float two = 2;
	    return r * ((((two / 9 * y + two / 7) * y + two / 5) * y +
			    two / 3) * y - x);
	} else {
	    static const Float tol_logcf = 1e-14;
	    return r * (2 * y * logcf ((Float)y, (Float)3, (Float)2, (Float)tol_logcf) - x);
	}
    }
}


/* Compute  log(gamma(a+1))  accurately also for small a (0 < a < 0.5). */
    template<class Float>
    Float lgamma1p (Float a)
{
    if (fabs (a) >= 0.5)
	return lgammafn (a + 1);

    const Float eulers_const =	 0.5772156649015328606065120900824024;

    /* coeffs[i] holds (zeta(i+2)-1)/(i+2) , i = 0:(N-1), N = 40 : */
    const int N = 40;
    static const Float coeffs[40] = {
	0.3224670334241132182362075833230126e-0,/* = (zeta(2)-1)/2 */
	0.6735230105319809513324605383715000e-1,/* = (zeta(3)-1)/3 */
	0.2058080842778454787900092413529198e-1,
	0.7385551028673985266273097291406834e-2,
	0.2890510330741523285752988298486755e-2,
	0.1192753911703260977113935692828109e-2,
	0.5096695247430424223356548135815582e-3,
	0.2231547584535793797614188036013401e-3,
	0.9945751278180853371459589003190170e-4,
	0.4492623673813314170020750240635786e-4,
	0.2050721277567069155316650397830591e-4,
	0.9439488275268395903987425104415055e-5,
	0.4374866789907487804181793223952411e-5,
	0.2039215753801366236781900709670839e-5,
	0.9551412130407419832857179772951265e-6,
	0.4492469198764566043294290331193655e-6,
	0.2120718480555466586923135901077628e-6,
	0.1004322482396809960872083050053344e-6,
	0.4769810169363980565760193417246730e-7,
	0.2271109460894316491031998116062124e-7,
	0.1083865921489695409107491757968159e-7,
	0.5183475041970046655121248647057669e-8,
	0.2483674543802478317185008663991718e-8,
	0.1192140140586091207442548202774640e-8,
	0.5731367241678862013330194857961011e-9,
	0.2759522885124233145178149692816341e-9,
	0.1330476437424448948149715720858008e-9,
	0.6422964563838100022082448087644648e-10,
	0.3104424774732227276239215783404066e-10,
	0.1502138408075414217093301048780668e-10,
	0.7275974480239079662504549924814047e-11,
	0.3527742476575915083615072228655483e-11,
	0.1711991790559617908601084114443031e-11,
	0.8315385841420284819798357793954418e-12,
	0.4042200525289440065536008957032895e-12,
	0.1966475631096616490411045679010286e-12,
	0.9573630387838555763782200936508615e-13,
	0.4664076026428374224576492565974577e-13,
	0.2273736960065972320633279596737272e-13,
	0.1109139947083452201658320007192334e-13/* = (zeta(40+1)-1)/(40+1) */
    };

    const Float c = 0.2273736845824652515226821577978691e-12;/* zeta(N+2)-1 */
    const Float tol_logcf = 1e-14;

    /* Abramowitz & Stegun 6.1.33 : for |x| < 2,
     * <==> log(gamma(1+x)) = -(log(1+x) - x) - gamma*x + x^2 * \sum_{n=0}^\infty c_n (-x)^n
     * where c_n := (Zeta(n+2) - 1)/(n+2)  = coeffs[n]
     *
     * Here, another convergence acceleration trick is used to compute
     * lgam(x) :=  sum_{n=0..Inf} c_n (-x)^n
     */
    Float lgam = c * logcf((Float)(-a / 2), (Float)(N + 2), (Float)1, (Float)tol_logcf);
    for (int i = N - 1; i >= 0; i--)
	lgam = coeffs[i] - a * lgam;

    return (a * lgam - eulers_const) * a - log1pmx (a);
} /* lgamma1p */



/*
 * Compute the log of a sum from logs of terms, i.e.,
 *
 *     log (exp (logx) + exp (logy))
 *
 * without causing overflows and without throwing away large handfuls
 * of accuracy.
 */
    template<class Float>
    Float logspace_add (Float logx, Float logy)
{
    return myfmax2 (logx, logy) + log1p (exp (-fabs (logx - logy)));
}


/*
 * Compute the log of a difference from logs of terms, i.e.,
 *
 *     log (exp (logx) - exp (logy))
 *
 * without causing overflows and without throwing away large handfuls
 * of accuracy.
 */
    template<class Float>
    Float logspace_sub (Float logx, Float logy)
{
    return logx + atomic::robust_utils::R_Log1_Exp(logy - logx);
}

/*
 * Compute the log of a sum from logs of terms, i.e.,
 *
 *     log (sum_i  exp (logx[i]) ) =
 *     log (e^M * sum_i  e^(logx[i] - M) ) =
 *     M + log( sum_i  e^(logx[i] - M)
 *
 * without causing overflows or throwing much accuracy.
 */
    template<class Float>
    Float logspace_sum (const Float* logx, int n)
{
    if(n == 0) return ML_NEGINF; // = log( sum(<empty>) )
    if(n == 1) return logx[0];
    if(n == 2) return logspace_add(logx[0], logx[1]);
    // else (n >= 3) :
    int i;
    // Mx := max_i log(x_i)
    Float Mx = logx[0];
    for(i = 1; i < n; i++) if(Mx < logx[i]) Mx = logx[i];
    Float s = (Float) 0.;
    for(i = 0; i < n; i++) s += exp(logx[i] - Mx);
    return Mx + (Float) log(s);
}



/* dpois_wrap (x__1, lambda) := dpois(x__1 - 1, lambda);  where
 * dpois(k, L) := exp(-L) L^k / gamma(k+1)  {the usual Poisson probabilities}
 *
 * and  dpois*(.., give_log = TRUE) :=  log( dpois*(..) )
*/
    template<class Float>
    static Float
dpois_wrap (Float x_plus_1, Float lambda, int give_log)
{
#ifdef DEBUG_p
    REprintf (" dpois_wrap(x+1=%.14g, lambda=%.14g, log=%d)\n",
	      x_plus_1, lambda, give_log);
#endif
    if (!R_FINITE(lambda))
	return R_D__0;
    if (x_plus_1 > 1)
      return dpois_raw ((Float)(x_plus_1 - 1), (Float)lambda, (int)give_log);
    if (lambda > fabs(x_plus_1 - 1) * M_cutoff)
	return R_D_exp(-lambda - lgammafn(x_plus_1));
    else {
      Float d = dpois_raw ((Float)x_plus_1, (Float)lambda, (int)give_log);
#ifdef DEBUG_p
	REprintf ("  -> d=dpois_raw(..)=%.14g\n", d);
#endif
	return give_log
	    ? d + log (x_plus_1 / lambda)
	    : d * (x_plus_1 / lambda);
    }
}

/*
 * Abramowitz and Stegun 6.5.29 [right]
 */
    template<class Float>
    static Float
pgamma_smallx (Float x, Float alph, int lower_tail, int log_p)
{
    Float sum = 0, c = alph, n = 0, term;

#ifdef DEBUG_p
    REprintf (" pg_smallx(x=%.12g, alph=%.12g): ", x, alph);
#endif

    /*
     * Relative to 6.5.29 all terms have been multiplied by alph
     * and the first, thus being 1, is omitted.
     */

    do {
	n = n+1;
	c *= -x / n;
	term = c / (alph + n);
	sum += term;
    } while (fabs (term) > DBL_EPSILON * fabs (sum));

#ifdef DEBUG_p
    REprintf ("%5.0f terms --> conv.sum=%g;", n, sum);
#endif
    if (lower_tail) {
	Float f1 = log_p ? log1p (sum) : 1 + sum;
	Float f2 = 0.0;
	if (alph > 1) {
	  f2 = dpois_raw ((Float)alph, (Float)x, (int)log_p);
	  f2 = log_p ? f2 + x : f2 * exp (x);
	} else if (log_p)
	    f2 = alph * log (x) - lgamma1p (alph);
	else
	    f2 = pow (x, alph) / exp (lgamma1p (alph));
#ifdef DEBUG_p
    REprintf (" (f1,f2)= (%g,%g)\n", f1,f2);
#endif
	return log_p ? f1 + f2 : f1 * f2;
    } else {
	Float lf2 = alph * log (x) - lgamma1p (alph);
#ifdef DEBUG_p
	REprintf (" 1:%.14g  2:%.14g\n", alph * log (x), lgamma1p (alph));
	REprintf (" sum=%.14g  log(1+sum)=%.14g	 lf2=%.14g\n",
		  sum, log1p (sum), lf2);
#endif
	if (log_p)
	    return atomic::robust_utils::R_Log1_Exp (log1p (sum) + lf2);
	else {
	    Float f1m1 = sum;
	    Float f2m1 = expm1 (lf2);
	    return -(f1m1 + f2m1 + f1m1 * f2m1);
	}
    }
} /* pgamma_smallx() */

    template<class Float>
    static Float
pd_upper_series (Float x, Float y, int log_p)
{
    Float term = x / y;
    Float sum = term;

    do {
	y = y+1;
	term *= x / y;
	sum += term;
    } while (term > sum * DBL_EPSILON);

    /* sum =  \sum_{n=1}^ oo  x^n     / (y*(y+1)*...*(y+n-1))
     *	   =  \sum_{n=0}^ oo  x^(n+1) / (y*(y+1)*...*(y+n))
     *	   =  x/y * (1 + \sum_{n=1}^oo	x^n / ((y+1)*...*(y+n)))
     *	   ~  x/y +  o(x/y)   {which happens when alph -> Inf}
     */
    return log_p ? log (sum) : sum;
}

/* Continued fraction for calculation of
 *    scaled upper-tail F_{gamma}
 *  ~=  (y / d) * [1 +  (1-y)/d +  O( ((1-y)/d)^2 ) ]
 */
    template<class Float>
    static Float
pd_lower_cf (Float y, Float d)
{
    Float f= 0.0 /* -Wall */, of, f0;
    Float i, c2, c3, c4,  a1, b1,  a2, b2;

#define	NEEDED_SCALE				\
	  (b2 > scalefactor) {			\
	    a1 /= scalefactor;			\
	    b1 /= scalefactor;			\
	    a2 /= scalefactor;			\
	    b2 /= scalefactor;			\
	}

#define max_it 200000

#ifdef DEBUG_p
    REprintf("pd_lower_cf(y=%.14g, d=%.14g)", y, d);
#endif
    if (y == 0) return 0;

    f0 = y/d;
    /* Needed, e.g. for  pgamma(10^c(100,295), shape= 1.1, log=TRUE): */
    if(fabs(y - 1) < fabs(d) * DBL_EPSILON) { /* includes y < d = Inf */
#ifdef DEBUG_p
	REprintf(" very small 'y' -> returning (y/d)\n");
#endif
	return (f0);
    }

    if(f0 > 1.) f0 = 1.;
    c2 = y;
    c4 = d; /* original (y,d), *not* potentially scaled ones!*/

    a1 = 0; b1 = 1;
    a2 = y; b2 = d;

    while NEEDED_SCALE

    i = 0; of = -1.; /* far away */
    while (i < max_it) {

	i=i+1;	c2=c2-1;	c3 = i * c2;	c4 += 2;
	/* c2 = y - i,  c3 = i(y - i),  c4 = d + 2i,  for i odd */
	a1 = c4 * a2 + c3 * a1;
	b1 = c4 * b2 + c3 * b1;

	i=i+1;	c2=c2-1;	c3 = i * c2;	c4 += 2;
	/* c2 = y - i,  c3 = i(y - i),  c4 = d + 2i,  for i even */
	a2 = c4 * a1 + c3 * a2;
	b2 = c4 * b1 + c3 * b2;

	if NEEDED_SCALE

	if (b2 != 0) {
	    f = a2 / b2;
	    /* convergence check: relative; "absolute" for very small f : */
	    if (fabs (f - of) <= DBL_EPSILON * myfmax2((Float)f0, (Float)fabs(f))) {
#ifdef DEBUG_p
		REprintf(" %g iter.\n", i);
#endif
		return f;
	    }
	    of = f;
	}
    }

    MATHLIB_WARNING(" ** NON-convergence in pgamma()'s pd_lower_cf() f= %g.\n",
		    f);
    return f;/* should not happen ... */
} /* pd_lower_cf() */
#undef NEEDED_SCALE


    template<class Float>
    static Float
pd_lower_series (Float lambda, Float y)
{
    Float term = 1, sum = 0;

#ifdef DEBUG_p
    REprintf("pd_lower_series(lam=%.14g, y=%.14g) ...", lambda, y);
#endif
    while (y >= 1 && term > sum * DBL_EPSILON) {
	term *= y / lambda;
	sum += term;
	y=y-1;
    }
    /* sum =  \sum_{n=0}^ oo  y*(y-1)*...*(y - n) / lambda^(n+1)
     *	   =  y/lambda * (1 + \sum_{n=1}^Inf  (y-1)*...*(y-n) / lambda^n)
     *	   ~  y/lambda + o(y/lambda)
     */
#ifdef DEBUG_p
    REprintf(" done: term=%g, sum=%g, y= %g\n", term, sum, y);
#endif

    if (y != floor (y)) {
	/*
	 * The series does not converge as the terms start getting
	 * bigger (besides flipping sign) for y < -lambda.
	 */
	Float f;
#ifdef DEBUG_p
	REprintf(" y not int: add another term ");
#endif
	/* FIXME: in quite few cases, adding  term*f  has no effect (f too small)
	 *	  and is unnecessary e.g. for pgamma(4e12, 121.1) */
	f = pd_lower_cf ((Float)y, (Float)(lambda + 1 - y));
#ifdef DEBUG_p
	REprintf("  (= %.14g) * term = %.14g to sum %g\n", f, term * f, sum);
#endif
	sum += term * f;
    }

    return sum;
} /* pd_lower_series() */

/*
 * Compute the following ratio with higher accuracy that would be had
 * from doing it directly.
 *
 *		 dnorm (x, 0, 1, FALSE)
 *	   ----------------------------------
 *	   pnorm (x, 0, 1, lower_tail, FALSE)
 *
 * Abramowitz & Stegun 26.2.12
 */
    template<class Float>
    static Float
dpnorm (Float x, int lower_tail, Float lp)
{
    /*
     * So as not to repeat a pnorm call, we expect
     *
     *	 lp == pnorm (x, 0, 1, lower_tail, TRUE)
     *
     * but use it only in the non-critical case where either x is small
     * or p==exp(lp) is close to 1.
     */

    if (x < 0) {
	x = -x;
	lower_tail = !lower_tail;
    }

    if (x > 10 && !lower_tail) {
	Float term = 1 / x;
	Float sum = term;
	Float x2 = x * x;
	Float i = 1;

	do {
	    term *= -i / x2;
	    sum += term;
	    i += 2;
	} while (fabs (term) > DBL_EPSILON * sum);

	return 1 / sum;
    } else {
      Float d = dnorm ((Float)x, (Float)0., (Float)1., (int)FALSE);
	return d / exp (lp);
    }
}

/*
 * Asymptotic expansion to calculate the probability that Poisson variate
 * has value <= x.
 * Various assertions about this are made (without proof) at
 * http://members.aol.com/iandjmsmith/PoissonApprox.htm
 */
    template<class Float>
    static Float
ppois_asymp (Float x, Float lambda, int lower_tail, int log_p)
{
    static const Float coefs_a[8] = {
	-1e99, /* placeholder used for 1-indexing */
	2/3.,
	-4/135.,
	8/2835.,
	16/8505.,
	-8992/12629925.,
	-334144/492567075.,
	698752/1477701225.
    };

    static const Float coefs_b[8] = {
	-1e99, /* placeholder */
	1/12.,
	1/288.,
	-139/51840.,
	-571/2488320.,
	163879/209018880.,
	5246819/75246796800.,
	-534703531/902961561600.
    };

    Float elfb, elfb_term;
    Float res12, res1_term, res1_ig, res2_term, res2_ig;
    Float dfm, pt_, s2pt, f, np;
    int i;

    dfm = lambda - x;
    /* If lambda is large, the distribution is highly concentrated
       about lambda.  So representation error in x or lambda can lead
       to arbitrarily large values of pt_ and hence divergence of the
       coefficients of this approximation.
    */
    pt_ = - log1pmx (dfm / x);
    s2pt = sqrt (2 * x * pt_);
    if (dfm < 0) s2pt = -s2pt;

    res12 = 0;
    res1_ig = res1_term = sqrt (x);
    res2_ig = res2_term = s2pt;
    for (i = 1; i < 8; i++) {
	res12 += res1_ig * coefs_a[i];
	res12 += res2_ig * coefs_b[i];
	res1_term *= pt_ / i ;
	res2_term *= 2 * pt_ / (2 * i + 1);
	res1_ig = res1_ig / x + res1_term;
	res2_ig = res2_ig / x + res2_term;
    }

    elfb = x;
    elfb_term = 1;
    for (i = 1; i < 8; i++) {
	elfb += elfb_term * coefs_b[i];
	elfb_term /= x;
    }
    if (!lower_tail) elfb = -elfb;
#ifdef DEBUG_p
    REprintf ("res12 = %.14g   elfb=%.14g\n", elfb, res12);
#endif

    f = res12 / elfb;

    np = pnorm5_raw ((Float)s2pt, (Float)0.0, (Float)1.0, (Float)!lower_tail, (Float)log_p);

    if (log_p) {
	Float n_d_over_p = dpnorm (s2pt, !lower_tail, np);
#ifdef DEBUG_p
	REprintf ("pp*_asymp(): f=%.14g	 np=e^%.14g  nd/np=%.14g  f*nd/np=%.14g\n",
		  f, np, n_d_over_p, f * n_d_over_p);
#endif
	return np + log1p (f * n_d_over_p);
    } else {
      Float nd = dnorm4 ((Float)s2pt, (Float)0., (Float)1., (int)log_p);

#ifdef DEBUG_p
	REprintf ("pp*_asymp(): f=%.14g	 np=%.14g  nd=%.14g  f*nd=%.14g\n",
		  f, np, nd, f * nd);
#endif
	return np + f * nd;
    }
} /* ppois_asymp() */


    template<class Float>
    Float pgamma_raw (Float x, Float alph, int lower_tail, int log_p)
{
/* Here, assume that  (x,alph) are not NA  &  alph > 0 . */

    Float res;

#ifdef DEBUG_p
    REprintf("pgamma_raw(x=%.14g, alph=%.14g, low=%d, log=%d)\n",
	     x, alph, lower_tail, log_p);
#endif
    R_P_bounds_01(x, 0., ML_POSINF);

    if (x < 1) {
	res = pgamma_smallx (x, alph, lower_tail, log_p);
    } else if (x <= alph - 1 && x < 0.8 * (alph + 50)) {
	/* incl. large alph compared to x */
	Float sum = pd_upper_series (x, alph, log_p);/* = x/alph + o(x/alph) */
	Float d = dpois_wrap ((Float)alph, (Float)x, (int)log_p);
#ifdef DEBUG_p
	REprintf(" alph 'large': sum=pd_upper*()= %.12g, d=dpois_w(*)= %.12g\n",
		 sum, d);
#endif
	if (!lower_tail)
	    res = log_p
		? atomic::robust_utils::R_Log1_Exp (d + sum)
		: 1 - d * sum;
	else
	    res = log_p ? sum + d : sum * d;
    } else if (alph - 1 < x && alph < 0.8 * (x + 50)) {
	/* incl. large x compared to alph */
	Float sum;
	Float d = dpois_wrap (alph, x, log_p);
#ifdef DEBUG_p
	REprintf(" x 'large': d=dpois_w(*)= %.14g ", d);
#endif
	if (alph < 1) {
	    if (x * DBL_EPSILON > 1 - alph)
		sum = R_D__1;
	    else {
	      Float f = pd_lower_cf ((Float)alph, (Float)(x - (alph - 1))) * x / alph;
		/* = [alph/(x - alph+1) + o(alph/(x-alph+1))] * x/alph = 1 + o(1) */
		sum = log_p ? log (f) : f;
	    }
	} else {
	  sum = pd_lower_series ((Float)x, (Float)(alph - 1));/* = (alph-1)/x + o((alph-1)/x) */
	    sum = log_p ? log1p (sum) : 1 + sum;
	}
#ifdef DEBUG_p
	REprintf(", sum= %.14g\n", sum);
#endif
	if (!lower_tail)
	    res = log_p ? sum + d : sum * d;
	else
	    res = log_p
		? atomic::robust_utils::R_Log1_Exp (d + sum)
		: 1 - d * sum;
    } else { /* x >= 1 and x fairly near alph. */
#ifdef DEBUG_p
	REprintf(" using ppois_asymp()\n");
#endif
	res = ppois_asymp ((Float)(alph - 1), (Float)x, (int)!lower_tail, (int)log_p);
    }

    /*
     * We lose a fair amount of accuracy to underflow in the cases
     * where the final result is very close to DBL_MIN.	 In those
     * cases, simply redo via log space.
     */
    if (!log_p && res < DBL_MIN / DBL_EPSILON) {
	/* with(.Machine, Float.xmin / double.eps) #|-> 1.002084e-292 */
#ifdef DEBUG_p
	REprintf(" very small res=%.14g; -> recompute via log\n", res);
#endif
	return exp (pgamma_raw (x, alph, lower_tail, 1));
    } else
	return res;
}

    template<class Float>
Float pgamma(Float x, Float alph, Float scale, int lower_tail, int log_p)
{
#ifdef IEEE_754
    if (ISNAN(x) || ISNAN(alph) || ISNAN(scale))
	return x + alph + scale;
#endif
    if(alph < 0. || scale <= 0.)
	ML_WARN_return_NAN;
    x /= scale;
#ifdef IEEE_754
    if (ISNAN(x)) /* eg. original x = scale = +Inf */
	return x;
#endif
    if(alph == 0.) /* limit case; useful e.g. in pnchisq() */
	return (x <= 0) ? R_DT_0: R_DT_1; /* <= assert  pgamma(0,0) ==> 0 */
    return pgamma_raw (x, alph, lower_tail, log_p);
}
/* From: terra@gnome.org (Morten Welinder)
 * To: R-bugs@biostat.ku.dk
 * Cc: maechler@stat.math.ethz.ch
 * Subject: Re: [Rd] pgamma discontinuity (PR#7307)
 * Date: Tue, 11 Jan 2005 13:57:26 -0500 (EST)

 * this version of pgamma appears to be quite good and certainly a vast
 * improvement over current R code.  (I last looked at 2.0.1)  Apart from
 * type naming, this is what I have been using for Gnumeric 1.4.1.

 * This could be included into R as-is, but you might want to benefit from
 * making logcf, log1pmx, lgamma1p, and possibly logspace_add/logspace_sub
 * available to other parts of R.

 * MM: I've not (yet?) taken  logcf(), but the other four
 */


