from biocomputedm.extensions import login_manager
from flask.ext.login import current_user

__author__ = "jon"

def login_required(*role):
    def wrapper(fn):
        from functools import wraps

        @wraps(fn)
        def decorated_view(*args, **kwargs):
            if current_user is not None and current_user.is_authenticated:
                if current_user.get_role() == "Site Admin":
                    return fn(*args, **kwargs)

                elif role == "ANY":
                    return fn(*args, **kwargs)

                else:
                    for item in role:
                        if item == "ANY":
                            return fn(*args, **kwargs)

                        if current_user.get_role() == item:
                            return fn(*args, **kwargs)

            return login_manager.unauthorized()

        return decorated_view

    return wrapper