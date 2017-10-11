""".. Ignore pydocstyle D400.

.. autoclass:: resolwe.test.ResolweAPITestCase
    :members:

"""

from __future__ import absolute_import, division, print_function, unicode_literals

import six

from django.contrib.auth import get_user_model
from django.contrib.contenttypes.models import ContentType
from django.core.urlresolvers import reverse
from django.test import override_settings

from rest_framework.test import APIRequestFactory, APITestCase, force_authenticate


@override_settings(CELERY_ALWAYS_EAGER=True)
class ResolweAPITestCase(APITestCase):
    """Base class for testing Resolwe REST API.

    This class is derived from Django REST Framework's
    :drf:`APITestCase <testing/#test-cases>` class and has implemented
    some basic features that makes testing Resolwe API easier. These
    features includes following functions:

    .. automethod:: _get_list
    .. automethod:: _get_detail
    .. automethod:: _post
    .. automethod:: _patch
    .. automethod:: _delete
    .. automethod:: _detail_permissions

    It also has included 2 views made from referenced DRF's ``ViewSet``.
    First mimic list view and has following links between request's
    methods and ViewSet's methods:

      *  ``GET`` -> ``list``
      *  ``POST`` -> ``create``

    Second mimic detail view and has following links between request's
    methods and ViewSet's methods:

      *  ``GET`` -> ``retrieve``
      *  ``PUT`` -> ``update``
      *  ``PATCH`` -> ``partial_update``
      *  ``DELETE`` -> ``destroy``
      *  ``POST`` -> ``permissions``

    If any of the listed methods is not defined in the VievSet,
    corresponding link is omitted.

    .. note::
        ``self.viewset`` (instance of DRF's ``Viewset``) and
        ``self.resource_name`` (string) must be defined before calling
        super ``setUp`` method to work properly.

    ``self.factory`` is instance of DRF's ``APIRequestFactory``.

    """

    def _pre_setup(self, *args, **kwargs):
        # NOTE: This is a work-around for Django issue #10827
        # (https://code.djangoproject.com/ticket/10827) that clears the
        # ContentType cache before permissions are setup.
        ContentType.objects.clear_cache()
        super(ResolweAPITestCase, self)._pre_setup(*args, **kwargs)

    def setUp(self):
        """Prepare data."""
        super(ResolweAPITestCase, self).setUp()

        # TODO: Remove this when removing fixtures
        if get_user_model().objects.filter(pk=2).exists():
            self.user1 = get_user_model().objects.get(pk=2)
        if get_user_model().objects.filter(pk=3).exists():
            self.user2 = get_user_model().objects.get(pk=3)
        if get_user_model().objects.filter(pk=4).exists():
            self.user3 = get_user_model().objects.get(pk=4)
        if get_user_model().objects.filter(pk=5).exists():
            self.admin = get_user_model().objects.get(pk=5)

        user_model = get_user_model()
        # TODO: Change username to `admin` when fixtures are removed
        self.admin = user_model.objects.create_superuser(username='admin2', email='admin@test.com', password='admin')
        self.contributor = user_model.objects.create_user(username='contributor')

        if not hasattr(self, 'viewset'):
            raise KeyError("`self.viewset` must be defined in child class")

        if not hasattr(self, 'resource_name'):
            raise KeyError("`self.resource_name` must be defined in child class")

        self.factory = APIRequestFactory()

        list_url_mapping = {}
        if hasattr(self.viewset, 'list'):  # pylint: disable=no-member
            list_url_mapping['get'] = 'list'
        if hasattr(self.viewset, 'create'):  # pylint: disable=no-member
            list_url_mapping['post'] = 'create'

        self.list_view = self.viewset.as_view(list_url_mapping)  # pylint: disable=no-member

        detail_url_mapping = {}
        if hasattr(self.viewset, 'retrieve'):  # pylint: disable=no-member
            detail_url_mapping['get'] = 'retrieve'
        if hasattr(self.viewset, 'update'):  # pylint: disable=no-member
            detail_url_mapping['put'] = 'update'
        if hasattr(self.viewset, 'partial_update'):  # pylint: disable=no-member
            detail_url_mapping['patch'] = 'partial_update'
        if hasattr(self.viewset, 'destroy'):  # pylint: disable=no-member
            detail_url_mapping['delete'] = 'destroy'
        if hasattr(self.viewset, 'detail_permissions'):  # pylint: disable=no-member
            detail_url_mapping['post'] = 'detail_permissions'

        self.detail_view = self.viewset.as_view(detail_url_mapping)  # pylint: disable=no-member

    def detail_url(self, pk):
        """Get detail url."""
        return reverse('resolwe-api:{}-detail'.format(self.resource_name), kwargs={'pk': pk})  # noqa pylint: disable=no-member

    def detail_permissions(self, pk):
        """Get detail permissions url."""
        return reverse('resolwe-api:{}-permissions'.format(self.resource_name), kwargs={'pk': pk})  # noqa pylint: disable=no-member

    @property
    def list_url(self):
        """Get list url."""
        return reverse('resolwe-api:{}-list'.format(self.resource_name))  # pylint: disable=no-member

    def _render_query_params(self, params):
        """Generate query parameters from given dict."""
        if not params:
            return ''

        return '?' + '&'.join('{}={}'.format(key, value) for key, value in six.iteritems(params))

    def _get_list(self, user=None, query_params={}):
        """Make ``GET`` request to ``self.list_view`` view.

        If ``user`` is not ``None``, the given user is authenticated
        before making the request.

        :param user: User to authenticate in request
        :type user: :class:`~django.contrib.auth.models.User` or :data:`None`
        :return: API response object
        :rtype: :drf:`Response <responses/#response>`

        """
        url = self.list_url + self._render_query_params(query_params)
        request = self.factory.get(url, format='json')
        force_authenticate(request, user)
        return self.list_view(request)

    def _get_detail(self, pk, user=None, query_params={}):
        """Make ``GET`` request to ``self.detail_view`` view.

        If ``user`` is not ``None``, the given user is authenticated
        before making the request.

        :param int pk: Primary key of the coresponding object
        :param user: User to authenticate in request
        :type user: :class:`~django.contrib.auth.models.User` or :data:`None`
        :return: API response object
        :rtype: :drf:`Response <responses/#response>`

        """
        url = self.detail_url(pk) + self._render_query_params(query_params)
        request = self.factory.get(url, format='json')
        force_authenticate(request, user)
        return self.detail_view(request, pk=pk)

    def _post(self, data={}, user=None, query_params={}):
        """Make ``POST`` request to ``self.list_view`` view.

        If ``user`` is not ``None``, the given user is authenticated
        before making the request.

        :param dict data: data for posting in request's body
        :param user: User to authenticate in request
        :type user: :class:`~django.contrib.auth.models.User` or :data:`None`
        :return: API response object
        :rtype: :drf:`Response <responses/#response>`

        """
        url = self.list_url + self._render_query_params(query_params)
        request = self.factory.post(url, data=data, format='json')
        force_authenticate(request, user)
        return self.list_view(request)

    def _patch(self, pk, data={}, user=None, query_params={}):
        """Make ``PATCH`` request to ``self.detail_view`` view.

        If ``user`` is not ``None``, the given user is authenticated
        before making the request.

        :param int pk: Primary key of the coresponding object
        :param dict data: data for posting in request's body
        :param user: User to authenticate in request
        :type user: :class:`~django.contrib.auth.models.User` or :data:`None`
        :return: API response object
        :rtype: :drf:`Response <responses/#response>`

        """
        url = self.detail_url(pk) + self._render_query_params(query_params)
        request = self.factory.patch(url, data=data, format='json')
        force_authenticate(request, user)
        return self.detail_view(request, pk=pk)

    def _delete(self, pk, user=None, query_params={}):
        """Make ``DELETE`` request to ``self.detail_view`` view.

        If ``user`` is not ``None``, the given user is authenticated
        before making the request.

        :param int pk: Primary key of the coresponding object
        :param user: User to authenticate in request
        :type user: :class:`~django.contrib.auth.models.User` or :data:`None`
        :return: API response object
        :rtype: :drf:`Response <responses/#response>`

        """
        url = self.detail_url(pk) + self._render_query_params(query_params)
        request = self.factory.delete(url, format='json')
        force_authenticate(request, user)
        return self.detail_view(request, pk=pk)

    def _detail_permissions(self, pk, data={}, user=None):
        """Make ``POST`` request to ``self.detail_view`` view.

        If ``user`` is not ``None``, the given user is authenticated
        before making the request.

        :param int pk: Primary key of the coresponding object
        :param dict data: data for posting in request's body
        :param user: User to authenticate in request
        :type user: :class:`~django.contrib.auth.models.User` or :data:`None`
        :return: API response object
        :rtype: :drf:`Response <responses/#response>`

        """
        request = self.factory.post(self.detail_permissions(pk), data=data, format='json')
        force_authenticate(request, user)
        return self.detail_view(request, pk=pk)

    def assertKeys(self, data, wanted):  # pylint: disable=invalid-name
        """Assert dictionary keys."""
        self.assertEqual(sorted(data.keys()), sorted(wanted))