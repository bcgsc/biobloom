/*
 * kseq_util.h
 *
 *  Created on: Jun 21, 2018
 *      Author: cjustin
 */

#ifndef KSEQ_UTIL_H_
#define KSEQ_UTIL_H_

#ifndef KSEQ_INIT_NEW
#define KSEQ_INIT_NEW
#include "kseq.h"
KSEQ_INIT(gzFile, gzread)
#endif /*KSEQ_INIT_NEW*/

static void cpy_kstr(kstring_t *dst, const kstring_t *src)
{
	if (src->l == 0) return;
	if (src->l + 1 > dst->m) {
		dst->m = src->l + 1;
		kroundup32(dst->m);
		dst->s = (char*)realloc(dst->s, dst->m);
	}
	dst->l = src->l;
	memcpy(dst->s, src->s, src->l + 1);
}


static void cpy_kseq(kseq_t *dst, const kseq_t *src)
{
	cpy_kstr(&dst->name, &src->name);
	cpy_kstr(&dst->seq,  &src->seq);
	cpy_kstr(&dst->qual, &src->qual);
	cpy_kstr(&dst->comment, &src->comment);
}

#endif /* KSEQ_UTIL_H_ */
