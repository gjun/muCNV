//
//  debug.h
//  muCNV
//
//  Created by Goo Jun on 11/13/18.
//  Copyright © 2018 Goo Jun. All rights reserved.
//

#ifndef debug_h
#define debug_h

#ifdef DEBUG
#define DMSG(str) do { std::cerr << str << std::endl; } while( false )
#else
#define DMSG(str) do { } while ( false )
#endif

#ifdef DDEBUG
#define DDMSG(str) do { std::cerr << str << std::endl; } while( false )
#define DDPRINT(fmt, ...) do { fprintf(stderr, fmt, __VA_ARGS__); } while ( false )
#else
#define DDMSG(str) do { } while ( false )
#define DDPRINT(fmt, ...) do { } while (false)
#endif

#endif /* debug_h */
