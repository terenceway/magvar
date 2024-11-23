/*
 * https://www.wayforward.net/magvar/av_geo.h
 * Copyright(c) 2005, 2008, 2022 Wayforward LLC
 *
 * Rather than having the geo-magnetic model in separate files, as
 * provided by NGDC, we code the data directly in C.
 *
 * ---
 *
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU Lesser General Public License as
 * published by the Free Software Foundation; either version 2 of the
 * License, or (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU Lesser General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public
 * License along with this program; if not, write to the Free Software
 * Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307
 * USA
 */
#ifndef AV_GEO_H
#define AV_GEO_H

#ifdef __cplusplus
extern "C" {
#endif

	struct GeoData {
		float x, y, z, t;
	};

	struct GeoModel {
		float yrmin, yrmax;
		float altmin, altmax;

		unsigned max1, max2, max3;

		const struct GeoData *data;
	};

	extern const struct GeoModel GEO_MODELS[];

	extern const unsigned GEO_MODEL_COUNT;

#ifdef __cplusplus
}
#endif

#endif /* AV_GEO_H */

/*
 * Overrides for Emacs so that we follow Linux's tabbing style.
 * Emacs will notice this stuff at the end of the file and automatically
 * adjust the settings for this buffer only.  This must remain at the
 * end of the file.
 * ----------------------------------------------------------------------
 * Local Variables:
 * c-file-style: "linux"
 * End:
 */
