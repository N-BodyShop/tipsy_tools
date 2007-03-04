/** @file xdrfile.hpp
 @brief routines for converting XDR format tipsy files.

 This code was originally part of ChaNGa and Salsa of the NChilada project.
 */

#include "xdr_template.hpp"

/// XDR conversion for the header structure
inline bool_t xdr_template(XDR* xdrs, tipsy_header* val) {
	return (xdr_template(xdrs, &(val->time))
			&& xdr_template(xdrs, &(val->nbodies))
			&& xdr_template(xdrs, &(val->ndim)) 
			&& xdr_template(xdrs, &(val->nsph))
			&& xdr_template(xdrs, &(val->ndark))
			&& xdr_template(xdrs, &(val->nstar))
			&& xdr_template(xdrs, &(val->zero)));
}

///XDR conversions for the particle types
inline bool_t xdr_template(XDR* xdrs, tipsy_gas_particle* p) {
	return (xdr_template(xdrs, &(p->mass))
		&& xdr_template(xdrs, &(p->pos[0]))
		&& xdr_template(xdrs, &(p->pos[1]))
		&& xdr_template(xdrs, &(p->pos[2]))
		&& xdr_template(xdrs, &(p->vel[0]))
		&& xdr_template(xdrs, &(p->vel[1]))
		&& xdr_template(xdrs, &(p->vel[2]))
		&& xdr_template(xdrs, &(p->rho))
		&& xdr_template(xdrs, &(p->temp))
		&& xdr_template(xdrs, &(p->hsmooth))
		&& xdr_template(xdrs, &(p->metals))
		&& xdr_template(xdrs, &(p->phi)));
}

inline bool_t xdr_template(XDR* xdrs, tipsy_dark_particle* p) {
	return (xdr_template(xdrs, &(p->mass))
		&& xdr_template(xdrs, &(p->pos[0]))
		&& xdr_template(xdrs, &(p->pos[1]))
		&& xdr_template(xdrs, &(p->pos[2]))
		&& xdr_template(xdrs, &(p->vel[0]))
		&& xdr_template(xdrs, &(p->vel[1]))
		&& xdr_template(xdrs, &(p->vel[2]))
		&& xdr_template(xdrs, &(p->eps))
		&& xdr_template(xdrs, &(p->phi)));
}

inline bool_t xdr_template(XDR* xdrs, tipsy_star_particle* p) {
	return (xdr_template(xdrs, &(p->mass))
		&& xdr_template(xdrs, &(p->pos[0]))
		&& xdr_template(xdrs, &(p->pos[1]))
		&& xdr_template(xdrs, &(p->pos[2]))
		&& xdr_template(xdrs, &(p->vel[0]))
		&& xdr_template(xdrs, &(p->vel[1]))
		&& xdr_template(xdrs, &(p->vel[2]))
		&& xdr_template(xdrs, &(p->metals))
		&& xdr_template(xdrs, &(p->tform))
		&& xdr_template(xdrs, &(p->eps))
		&& xdr_template(xdrs, &(p->phi)));
}
