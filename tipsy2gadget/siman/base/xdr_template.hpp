/** @file xdr_template.h
 Provides inline functions xdr_template() which perform
 XDR conversion of a value.  Numerous overloads
 of this function are provided for many common types.
 If you attempt to use this function for a type that does not have
 a provided overload, you will get a compile time error.
 You can provide overloads for your own types by copying the
 syntax used below.
 This framework allows you to call xdr_template() whenever you wish
 to do XDR conversion, and not worry about the underlying calls
 to xdr_int, xdr_float, etc.
 Originally written by Graeme Lufkin as part of ChaNGa and Salsa
 (released under GPL)
 */
#ifndef XDR_TEMPLATE_HPP
#define XDR_TEMPLATE_HPP

#include <rpc/rpc.h>

inline bool_t xdr_template(XDR* xdrs, int* val) {
	return xdr_int(xdrs, val);
}

inline bool_t xdr_template(XDR* xdrs, float* val) {
	return xdr_float(xdrs, val);
}

inline bool_t xdr_template(XDR* xdrs, double* val) {
	return xdr_double(xdrs, val);
}

#endif
