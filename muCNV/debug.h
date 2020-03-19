//
//  debug.h
//  muCNV
//
//  Created by Goo Jun on 11/13/18.
//  Copyright © 2018 Goo Jun. All rights reserved.
//

#ifndef debug_h
#define debug_h

#ifdef DDEBUG
#define DDEBUG_TEST 1
#else
#define DDEBUG_TEST 0
#endif

#ifdef DEBUG
#define DMSG(str) do { std::cerr << str << std::endl; } while( false )
#else
#define DMSG(str) do { } while ( false )
#endif

#define DDMSG(str) do { if (DDEBUG_TEST) std::cerr << str << std::endl; } while( false )
#define DDPRINT(fmt, ...) do { if (DDEBUG_TEST) fprintf(stderr, fmt, __VA_ARGS__); } while ( false )


#endif /* debug_h */
