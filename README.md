# Lattice-Aided Extraction of Spread-Spectrum Hidden Data

**Authors:** Fan Yang, Shanxiang Lyu, Hao Cheng, Jinming Wen, Hao Chen
*College of Cyber Security/College of Information Science and Technology, Jinan University, Guangzhou 510632, China*
*State Key Laboratory of Cryptology, P. O. Box 5159, Beijing, 100878, China*

## Abstract
This paper delves into the challenges of spread spectrum (SS) watermarking extraction, considering both blind and non-blind extraction scenarios, within the framework of lattice decoding. The orthogonality of carriers plays a crucial role in the accuracy of extraction, impacting the bit error rate (BER). When carriers lack sufficient orthogonality, conventional blind extraction methods such as multi-carrier iterative generalized least-squares (M-IGLS) and non-blind extraction techniques like MMSE-based schemes encounter performance degradation, posing difficulties in accurately recovering hidden data at the receiver end. To address these challenges, we propose two novel SS watermarking extraction approaches by integrating precise lattice decoding algorithms. Firstly, we introduce the highly accurate yet computationally efficient successive interference cancellation (SIC) algorithm to augment M-IGLS, resulting in a new method termed multi-carrier iterative successive interference cancellation (M-ISIC). Secondly, we adapt the near-optimal sphere decoding (SD) technique for non-blind extraction in SS watermarking. Theoretical analysis and experimental simulations showcase that our proposed M-ISIC and SD methods outperform M-IGLS and MMSE-based detectors, particularly in scenarios where carrier orthogonality is limited, achieving lower BER.

## MATLAB Codes
The MATLAB codes compare different extracting algorithms including:

- MMSE (Minimum Mean Square Error)
- JADE (Jean-Franccois Cardoso, et al. "Blind beamforming for non Gaussian signals", "High-order contrasts for independent component analysis")
- M_IGLS (Iterative General Least Square)
- M_ISIC (Iterative Successive Interference Cancellation)
- sph_dec (Sphere Decoding)

The main file of blind extraction is: `blind_extra_main.m`.
The main file of non-blind extraction is: `non_blind_extra_main.m`.
