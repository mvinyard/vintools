def _report_imports(env_globals=globals(), return_list=False, silent=False, remove_native_modules = ['__builtin__', '__builtins__', 'types', 'sys']):
    
    """
    Parameters:
    -----------
    return_list [ optional ]
        type:bool
    
    silent [ optional ]
        type:bool
    
    remove_native_modules [ optional ]
        type:list
    
    Returns:
    --------
    imported_modules [ optional ]
        type: list
    """
    
    from types import ModuleType
    
    imported_modules = []
    for name, val in env_globals.items():
        if name not in remove_native_modules:
            if isinstance(val, ModuleType):
                imported_modules.append(name)
                if not silent:
                    print(name)
    if return_list:
        return imported_modules