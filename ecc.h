#ifndef _ECC_H_
#define _ECC_H_

/* 
 * Enhanced Correlation Coefficient (ECC) Image Alignment Algorithm 
 * http://xanthippi.ceid.upatras.gr/people/evangelidis/ecc
 *
 * Authors: Georgios Evangelidis
 * e-mail: evagelid@ceid.upatras.gr
 * e-mail: georgios.evangelidis@inria.fr
 *
 * Copyright (2010): G.D.Evangelidis
 *
 * For details, see the paper:
 * G.D. Evangelidis, E.Z. Psarakis, "Parametric Image Alignment using 
 * Enhanced Correlation Coefficient Maximization, IEEE Transaction on
 * Pattern Analysis & Machine Intelligence, vol. 30, no. 10, 2008
 */

#include<opencv\cxcore.h>

enum WARP_MODE
{ 
	WARP_MODE_TRANSLATION,
	WARP_MODE_EUCLIDEAN,
	WARP_MODE_AFFINE,
	WARP_MODE_HOMOGRAPHY
};

/* This might be a suitable help text:
 
 findTransform
 ----------------------
 Calculates a translation, euclidean, affine transform or homography from a pair of images.
 
 .. ocv:function:: Mat findTransform( const Mat src_image, const Mat dst_image, const WARP_MODE mode, Mat map_matrix )
 
 .. ocv:pyfunction:: cv2.findTransform(src_image, dst_image, const WARP_MODE mode, map_matrix) -> retval
 
 .. ocv:cfunction:: CvMat* cvFindTransform( const IplImage* src_image, const IplImage* dst_image, const WARP_MODE mode, CvMat* map_matrix )
  
 :param src_image: Source image that is projected onto destination image.
 
 :param dst_image: Destination image.
 
 :param mode: pass TRANSLATION, EUCLIDEAN, AFFINE or HOMOGRAPHY to select mapping type
 
 :param map_matrix: Matrix holding initial transform.
 
 The function updates the :math:`2 \times 3` matrix of a {affine, euclidean, affine} transform or the :math:`3 \times 3` matrix of a homography so that
 
 .. math::
 
 \begin{bmatrix} x'_i \\ y'_i \end{bmatrix} = \texttt{map\_matrix} \cdot \begin{bmatrix} x_i \\ y_i \\ 1 \end{bmatrix}
 
 where
 
 .. math::
 
 src_image(i)=(x_i, y_i),
 dst_image(i)=(x'_i,y'_i),
 i=0,1,2,...,k
 
 .. seealso::
 
 :ocv:func:`warpAffine`,
 :ocv:func:`transform`
 */
CvMat * cvFindTransform(const IplImage * src_image, const IplImage * dst_image, CvMat * map_matrix, const WARP_MODE warp_mode = WARP_MODE_AFFINE, const CvTermCriteria & termination = cvTermCriteria (CV_TERMCRIT_ITER+CV_TERMCRIT_EPS, 50, 0.001));

#endif //_ECC_H_
