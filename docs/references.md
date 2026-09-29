# References

Access date for every online resource: **2026-09-29**.

Verification legend (how *we* checked each entry; nothing here is cited from memory alone):

* **FT** – full text (or the cited sections) obtained and read; section/equation numbers quoted in
  `theory_review.md` were read off the page images of the PDF.
* **FT-partial** – full text obtained; only the sections named were read.
* **ABS** – bibliographic metadata verified via the Crossref REST API
  (`https://api.crossref.org/works/<doi>`) **and** the abstract read (via OpenAlex
  `https://api.openalex.org/works/doi:<doi>` or Semantic Scholar). Claims attributed to these
  entries are restricted to what the abstract states.
* **META** – bibliographic metadata verified only (Crossref / publisher / catalogue); content not
  read. Such entries are cited only for existence/attribution, never for technical detail.

Where a claim in `theory_review.md` is our own derivation, it is labelled *(project analysis)*
and is not attributed to any of these works.

---

## Primary NRW sources

**[R1]** A. M. Nicolson and G. F. Ross, "Measurement of the intrinsic properties of materials by
time-domain techniques," *IEEE Trans. Instrum. Meas.*, vol. IM-19, no. 4, pp. 377–382, Nov. 1970.
DOI: 10.1109/TIM.1970.4313932. URL: https://doi.org/10.1109/TIM.1970.4313932.
*Verification:* ABS (Crossref + OpenAlex abstract). Full text not accessed (paywalled). The NRW
equations used in this project are taken from their reproduction in [R4], [R6], [R7].

**[R2]** W. B. Weir, "Automatic measurement of complex dielectric constant and permeability at
microwave frequencies," *Proc. IEEE*, vol. 62, no. 1, pp. 33–36, Jan. 1974.
DOI: 10.1109/PROC.1974.9382. URL: https://doi.org/10.1109/PROC.1974.9382.
*Verification:* ABS. Full text not accessed. The group-delay branch-selection procedure is
attributed to Weir by [R4] §2.1 and [R6] §2.2.1; we quote it from those reproductions.

## Baker-Jarvis / NIST formulations

**[R3]** J. Baker-Jarvis, E. J. Vanzura, and W. A. Kissick, "Improved technique for determining
complex permittivity with the transmission/reflection method," *IEEE Trans. Microw. Theory Tech.*,
vol. 38, no. 8, pp. 1096–1103, Aug. 1990. DOI: 10.1109/22.57336.
URL: https://doi.org/10.1109/22.57336.
*Verification:* ABS. **Full text NOT read** (IEEE paywall; the Semantic Scholar open-access link
returned HTTP 202/empty on repeated attempts). The abstract states: "The special case of
permittivity measurement is examined in detail. Robust algorithms for permittivity determination
that eliminate the ill-behaved nature of the commonly used procedures at frequencies corresponding
to integer multiples of one-half wavelength in the sample are presented." Equation-level content is
taken from the companion NIST report [R4] by the same first author in the same year, whose abstract
is nearly identical; the equation numbering of the journal paper is therefore **not** given.

**[R4]** J. Baker-Jarvis, "Transmission/reflection and short-circuit line permittivity
measurements," NIST Technical Note 1341, National Institute of Standards and Technology, Boulder,
CO, July 1990, 154 pp. DOI: 10.6028/NIST.TN.1341.
URL: https://nvlpubs.nist.gov/nistpubs/Legacy/TN/nbstechnicalnote1341.pdf.
*Verification:* FT-partial (Crossref DOI resolves to this PDF; read: Abstract, Ch. 1, §2.1
eqs. (2.12)–(2.37), §2.2 eq. (2.38) and text, Appendix E (Newton–Raphson, eqs. (6.177)–(6.180)),
Appendix F (Kramers–Kronig, eqs. (6.181)–(6.184))).

**[R5]** J. Baker-Jarvis, M. D. Janezic, J. H. Grosvenor Jr., and R. G. Geyer,
"Transmission/reflection and short-circuit line methods for measuring permittivity and
permeability," NIST Technical Note 1355, May 1992, 120 pp. DOI: 10.6028/NIST.TN.1355.
URL: https://nvlpubs.nist.gov/nistpubs/Legacy/TN/nbstechnicalnote1355.pdf.
*Verification:* META + PDF obtained and table of contents / §2.2 compared by text search with
[R6] (same section structure and same NRW and reference-plane-invariant equations). We cite the
revised edition [R6] for all technical content.

**[R6]** J. Baker-Jarvis, M. D. Janezic, J. H. Grosvenor Jr., and R. G. Geyer,
"Transmission/reflection and short-circuit line methods for measuring permittivity and
permeability," NIST Technical Note 1355-R (revised), Dec. 1993, 236 pp. (as printed on the title
verso). DOI: 10.6028/NIST.TN.1355r.
URL: https://nvlpubs.nist.gov/nistpubs/Legacy/TN/nbstechnicalnote1355r.pdf.
*Verification:* FT-partial. Read: Abstract; Ch. 1; §2.1.3 eqs. (2.29)–(2.39); §2.2.1–2.2.3
eqs. (2.40)–(2.61); §2.3 (intro); §2.5 eqs. (2.67)–(2.68) and surrounding text; Ch. 3 (§3.1–3.5,
optimized/regression solution); Ch. 5 (Discussion, Tables 5.1 and 5.2); Appendix C
(gap corrections, eqs. (C.16)–(C.34), §C.5); Appendix D.2–D.3 (eqs. (D.6), (D.8), (D.9)).

**[R7]** J. Baker-Jarvis, M. D. Janezic, B. F. Riddle, R. T. Johnk, P. Kabos, C. L. Holloway,
R. G. Geyer, and C. A. Grosvenor, "Measuring the permittivity and permeability of lossy materials:
solids, liquids, metals, building materials, and negative-index materials," NIST Technical Note
1536, Feb. 2005. DOI: 10.6028/NIST.TN.1536.
URL: https://nvlpubs.nist.gov/nistpubs/Legacy/TN/nbstechnicalnote1536.pdf
(also https://www.nist.gov/publications/measuring-permittivity-and-permeability-lossy-materials-solids-liquids-metals-and).
*Verification:* FT-partial (read §7.3.1–7.3.2 eqs. (97)–(101), §7.5.1–7.5.5 eqs. (112)–(128) and
Fig. 26). Author list checked against the PDF title page (includes P. Kabos; the NIST
publication web page omits him). 160-page scanned PDF.

**[R8]** J. Baker-Jarvis, R. G. Geyer, and P. D. Domich, "A nonlinear least-squares solution with
causality constraints applied to transmission line permittivity and permeability determination,"
*IEEE Trans. Instrum. Meas.*, vol. 41, no. 5, pp. 646–652, 1992. DOI: 10.1109/19.177336.
*Verification:* ABS. (Technical content cited from its restatement in [R6] Ch. 3, which was read.)

**[R9]** P. D. Domich, J. Baker-Jarvis, and R. G. Geyer, "Optimization techniques for permittivity
and permeability determination," NIST Interagency Report NISTIR 4571, 1991.
DOI: 10.6028/NIST.IR.4571. *Verification:* META (Crossref). Cited as ref. [36] in [R6].

**[R10]** A.-H. Boughriet, C. Legrand, and A. Chapoton, "Noniterative stable transmission/reflection
method for low-loss material complex permittivity determination," *IEEE Trans. Microw. Theory
Tech.*, vol. 45, no. 1, pp. 52–57, Jan. 1997. DOI: 10.1109/22.552032.
*Verification:* ABS. **Full text NOT read** (paywalled, no open copy found). The formula attributed
to this paper in `theory_review.md` is taken from NIST's restatement in [R7] eq. (99), which
explicitly names Boughriet; the ε_eff/μ_eff notation is additionally cross-checked against the
open-source implementation [R36] (secondary). The abstract states the method is for permittivity,
is "based on a simplified version of the well-known Nicolson-Ross-Weir (NRW) method", and shows "no
divergence ... at frequencies corresponding to integer multiples of one half wavelength" for
low-loss materials.

**[R11]** S. Kim and J. Baker-Jarvis, "An approximate approach to determining the permittivity and
permeability near λ/2 resonances in transmission/reflection measurements," *Progress In
Electromagnetics Research B*, vol. 58, pp. 95–109, 2014. DOI: 10.2528/PIERB13121308.
URL (NIST copy): https://tsapps.nist.gov/publication/get_pdf.cfm?pub_id=914263.
*Verification:* FT-partial (Abstract, §1, §2.1 incl. eqs. (1)–(16) and discussion of the
"BJ iterative method" eqs. (12)–(13)).

**[R12]** K. Chalapat, K. Sarvala, J. Li, and G. S. Paraoanu, "Wideband reference-plane invariant
method for measuring electromagnetic parameters of materials," *IEEE Trans. Microw. Theory Tech.*,
vol. 57, no. 9, pp. 2257–2267, Sep. 2009. DOI: 10.1109/TMTT.2009.2027160.
Author-posted reprint: arXiv:1009.4277, https://arxiv.org/abs/1009.4277.
*Verification:* FT-partial (arXiv reprint; §I–III read). Crossref confirms bibliographic data.

**[R31]** E. J. Vanzura, J. Baker-Jarvis, J. H. Grosvenor, and M. D. Janezic, "Intercomparison of
permittivity measurements using the transmission/reflection method in 7-mm coaxial transmission
lines," *IEEE Trans. Microw. Theory Tech.*, vol. 42, no. 11, pp. 2063–2070, 1994.
DOI: 10.1109/22.330120. *Verification:* ABS.

## Phase-branch / retrieval literature

**[R13]** U. C. Hasar and C. R. Westgate, "A broadband and stable method for unique complex
permittivity determination of low-loss materials," *IEEE Trans. Microw. Theory Tech.*, vol. 57,
no. 2, pp. 471–477, Feb. 2009. DOI: 10.1109/TMTT.2008.2011242. *Verification:* ABS.

**[R14]** J. J. Barroso and U. C. Hasar, "Resolving phase ambiguity in the inverse problem of
transmission/reflection measurement methods," *J. Infrared Millim. Terahertz Waves*, vol. 32,
no. 6, pp. 857–866, 2011. DOI: 10.1007/s10762-011-9792-7. *Verification:* META only (no abstract
available to us); cited for existence only.

**[R15]** V. V. Varadan and R. Ro, "Unique retrieval of complex permittivity and permeability of
dispersive materials from reflection and transmitted fields by enforcing causality," *IEEE Trans.
Microw. Theory Tech.*, vol. 55, no. 10, pp. 2224–2230, Oct. 2007. DOI: 10.1109/TMTT.2007.906473.
*Verification:* ABS.

**[R16]** Z. Szabó, G.-H. Park, R. Hedge, and E.-P. Li, "A unique extraction of metamaterial
parameters based on Kramers–Kronig relationship," *IEEE Trans. Microw. Theory Tech.*, vol. 58,
no. 10, pp. 2646–2653, Oct. 2010. DOI: 10.1109/TMTT.2010.2065310. *Verification:* ABS.
*Note:* this paper is in IEEE Trans. MTT, **not** Phys. Rev. B as sometimes mis-cited.

**[R17]** S. Arslanagić, T. V. Hansen, N. A. Mortensen, A. H. Gregersen, O. Sigmund,
R. W. Ziolkowski, and O. Breinbjerg, "A review of the scattering-parameter extraction method with
clarification of ambiguity issues in relation to metamaterial homogenization," *IEEE Antennas
Propag. Mag.*, vol. 55, no. 2, pp. 91–106, Apr. 2013. DOI: 10.1109/MAP.2013.6529320.
Preprint: arXiv:1206.4801. *Verification:* ABS + FT-partial (arXiv §5 "branch ambiguity" read).

**[R18]** D. R. Smith, S. Schultz, P. Markoš, and C. M. Soukoulis, "Determination of effective
permittivity and permeability of metamaterials from reflection and transmission coefficients,"
*Phys. Rev. B*, vol. 65, 195104, 2002. DOI: 10.1103/PhysRevB.65.195104. *Verification:* ABS.

**[R19]** D. R. Smith, D. C. Vier, T. Koschny, and C. M. Soukoulis, "Electromagnetic parameter
retrieval from inhomogeneous metamaterials," *Phys. Rev. E*, vol. 71, 036617, 2005.
DOI: 10.1103/PhysRevE.71.036617. *Verification:* ABS.

**[R20]** X. Chen, T. M. Grzegorczyk, B.-I. Wu, J. Pacheco Jr., and J. A. Kong, "Robust method to
retrieve the constitutive effective parameters of metamaterials," *Phys. Rev. E*, vol. 70, 016608,
2004. DOI: 10.1103/PhysRevE.70.016608. *Verification:* ABS.

**[R21]** T. Koschny, P. Markoš, D. R. Smith, and C. M. Soukoulis, "Resonant and antiresonant
frequency dependence of the effective parameters of metamaterials," *Phys. Rev. E*, vol. 68,
065602, 2003. DOI: 10.1103/PhysRevE.68.065602. *Verification:* ABS.

**[R22]** O. Luukkonen, S. I. Maslovski, and S. A. Tretyakov, "A stepwise Nicolson–Ross–Weir-based
material parameter extraction method," *IEEE Antennas Wireless Propag. Lett.*, vol. 10,
pp. 1295–1298, 2011. DOI: 10.1109/LAWP.2011.2175897. *Verification:* ABS.

**[R23]** E. J. Rothwell, J. L. Frasch, S. M. Ellison, P. Chahal, and R. O. Ouedraogo, "Analysis of
the Nicolson-Ross-Weir method for characterizing the electromagnetic properties of engineered
materials," *Progress In Electromagnetics Research*, vol. 157, pp. 31–47, 2016.
DOI: 10.2528/PIER16071706. *Verification:* ABS.

**[R30]** L. F. Chen, C. K. Ong, C. P. Neo, V. V. Varadan, and V. K. Varadan, *Microwave
Electronics: Measurement and Materials Characterization*. Chichester, U.K.: Wiley, 2004.
DOI: 10.1002/0470020466. *Verification:* META (Crossref).

## Causality, Kramers–Kronig, passivity

**[R24]** G. W. Milton, D. J. Eyre, and J. V. Mantese, "Finite frequency range Kramers-Kronig
relations: bounds on the dispersion," *Phys. Rev. Lett.*, vol. 79, no. 16, pp. 3062–3065, 1997.
DOI: 10.1103/PhysRevLett.79.3062. *Verification:* ABS.

**[R25]** V. Lucarini, J. J. Saarinen, K.-E. Peiponen, and E. M. Vartiainen, *Kramers-Kronig
Relations in Optical Materials Research*, Springer Series in Optical Sciences, vol. 110. Berlin:
Springer, 2005. DOI: 10.1007/b138913. ISBN 3-540-23673-2.
*Verification:* META (Crossref DOI/ISBN; author list from publisher/retailer listings). Content
not read; cited only as a general monograph.

**[R26]** K. F. Palmer, M. Z. Williams, and B. A. Budde, "Multiply subtractive Kramers–Kronig
analysis of optical data," *Appl. Opt.*, vol. 37, no. 13, p. 2660 ff., 1998.
DOI: 10.1364/AO.37.002660. *Verification:* ABS (Crossref gives the first page only; last page not verified).

**[R27]** R. Z. Bachrach and F. C. Brown, "Exciton-optical properties of TlBr and TlCl," *Phys.
Rev. B*, vol. 1, no. 2, pp. 818–831, 1970. DOI: 10.1103/PhysRevB.1.818.
*Verification:* ABS. *Caveat:* the abstract does not mention Kramers–Kronig analysis. The
attribution of the singly-subtractive KK relation to this paper follows common citation practice
(e.g. the MSKK literature [R26]) and was **not** verified from its full text.

**[R28]** P. Triverio, S. Grivet-Talocia, M. S. Nakhla, F. G. Canavero, and R. Achar, "Stability,
causality, and passivity in electrical interconnect models," *IEEE Trans. Adv. Packag.*, vol. 30,
no. 4, pp. 795–808, Nov. 2007. DOI: 10.1109/TADVP.2007.901567. *Verification:* ABS.

**[R29]** L. D. Landau, E. M. Lifshitz, and L. P. Pitaevskii, *Electrodynamics of Continuous
Media*, 2nd ed. (Course of Theoretical Physics, vol. 8). Oxford: Pergamon, 1984.
ISBN 978-0-08-030275-1. Chapter IX "The electromagnetic wave equations", pp. 257–289,
DOI: 10.1016/B978-0-08-030275-1.50015-1. *Verification:* META (Crossref chapter record). Content
not re-read for this project; the individual § numbers for the dispersion/analyticity/μ(ω)
discussion are **not** verified and are therefore not quoted.

## Uncertainty

**[R32]** JCGM 100:2008, *Evaluation of measurement data — Guide to the expression of uncertainty
in measurement* (GUM 1995 with minor corrections). BIPM, 2008.
URL: https://www.bipm.org/documents/20126/2071204/JCGM_100_2008_E.pdf.
*Verification:* PDF retrieved, title page read.

**[R33]** JCGM 101:2008, *Evaluation of measurement data — Supplement 1 to the GUM — Propagation of
distributions using a Monte Carlo method*. BIPM, 2008.
URL: https://www.bipm.org/documents/20126/2071204/JCGM_101_2008_E.pdf.
*Verification:* PDF retrieved, title page read.

## Instrument-vendor documentation

**[R34]** Keysight Technologies, "Keysight 85071E Materials Measurement Software — Technical
Overview," literature no. 5988-9472EN.
URL: https://www.keysight.com/us/en/assets/7018-01155/technical-overviews/5988-9472.pdf.
*Verification:* FT (algorithm table read). Lists "Reflection/Transmission Mu and Epsilon
(Nicholson-Ross-Weir, NRW)" → εr, μr; "Reflection/Transmission Epsilon Precision (NIST
Precision)" → εr only; "Transmission Epsilon Fast"; polynomial-fit (Bartley–Begley) models;
short-/arbitrary-backed reflection models.

**[R35]** Keysight Technologies, "N1500A Materials Measurement Suite — Technical Overview,"
literature no. 5992-0263EN.
URL: https://www.keysight.com/content/dam/keysight/en/doc/ungate/technical-overviews/5992-0263.pdf.
*Verification:* FT (algorithm table read). States that "NIST Precision" returns εr only, is "best
for longer samples of low-loss dielectric materials", and cites [R3] as its reference; states the
two μr = 1 two-port algorithms "do not suffer from discontinuities at frequencies where the sample
length is a multiple of half-wavelengths like the Nicholson-Ross model".

**[R37]** Agilent/Keysight Technologies, "Basics of Measuring the Dielectric Properties of
Materials," Application Note 5989-2589EN. Versions read: Agilent, printed June 26, 2006
(copy at https://anlage.umd.edu/5989-2589EN.pdf) and Keysight, published Aug. 20, 2020
(copy at https://mriquestions.com/uploads/3/4/5/7/34572113/5989-2589.pdf; official asset page
https://www.keysight.com/us/en/assets/7018-01284/application-notes/5989-2589.pdf).
*Verification:* FT (text search). Both versions mention the "precision NIST algorithm" in the
waveguide example; **neither version contains the phrase "New non-iterative" or "NIST
Precision (iterative)"** — the algorithm names asked about were not found in this app note (they
may exist in the 85071 on-line help, which we could not access). This is recorded as an
unverified attribution.

## Secondary / not independently verified

**[R36]** A. Boivin, *permittivitycalc* (open-source Python package), GitHub repository
https://github.com/boivinalex/permittivitycalc, file `permittivitycalc/sparam_data.py`,
method `_permittivity_calc`. *Verification:* source code read. Secondary source only: it implements
the "New Non-iterative method" of [R10] as ε_eff = (λ_og/Λ)², μ_eff = 1,
ε_r = (1 − λ0²/λc²) ε_eff + (λ0²/λc²)/μ_eff. (The repository line
`(lam_0**2/LAM_C*2)/mu_eff` contains an evident typo, `*2` for `**2`; irrelevant for coax where
λc → ∞.)

**[R38]** W. P. Westphal, "Techniques of measuring the permittivity and permeability of liquids and
solids in the frequency range 3 c/s to 50 kmc/s," Laboratory for Insulation Research Technical
Report, MIT, 1950. *Verification:* **NOT verified independently**; known to us only as ref. [34]
of [R6]. Cited solely as "gap-correction work cited by [R6]".

**[R39]** J. Baker-Jarvis, M. D. Janezic, and J. Krupka, "Measurements of coaxial dielectric
samples employing both transmission/reflection and resonant techniques to enhance air-gap
corrections," *2006 Int. Conf. Microwaves, Radar & Wireless Communications (MIKON)*,
pp. 1093–1096, 2006. DOI: 10.1109/MIKON.2006.4345376. *Verification:* META (Crossref).

---

### Items requested but not verified / corrected

* "Szabó et al. 2010 PRB" — the verified paper is [R16] in IEEE Trans. MTT.
* The phrases "NIST Precision (iterative)" and "New non-iterative" as 85071E algorithm names: the
  85071E technical overview [R34] lists "NIST Precision" (εr only) but we did **not** find
  "New non-iterative" in [R34], [R35] or either version of [R37]. A web search snippet suggests
  the 85071 on-line help cites Boughriet [R10], but we could not open that help file.
* Weir's and Nicolson–Ross's original full texts, the Baker-Jarvis 1990 journal paper and the
  Boughriet 1997 paper were not accessible; all equations attributed to them are quoted from the
  NIST reproductions [R4], [R6], [R7] and labelled as such.
