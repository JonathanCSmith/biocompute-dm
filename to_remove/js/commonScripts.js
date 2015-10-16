var checkAdminStatus = function() {
    var cookies = document.cookie;
    var lgn = cookies.indexOf("Exologin");

    if (lgn != -1) {
        $("#login").hide();
        $("#logout").show();
    }

    else {
        $("#login").show();
        $("#logout").hide();
    }
};

var main = function() {
    checkAdminStatus();
};

$(document).ready(main());